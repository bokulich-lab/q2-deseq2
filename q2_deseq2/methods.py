# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
from pathlib import Path
from subprocess import CalledProcessError, run
import tempfile

import biom
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype

from q2_deseq2._run_data import DESeq2RunResult, write_run_result_artifact
from q2_deseq2.types import DESeq2RunDirectoryFormat

_FORMULA_ALLOWED_CHARS = re.compile(r"^[A-Za-z0-9_:+*() \t]+$")
_FORMULA_TOKEN_RE = re.compile(r"\b[A-Za-z_][A-Za-z0-9_]*\b")
_COEF_SPEC_RE = re.compile(r"^coef::(.+)$")
_CONTRAST_SPEC_RE = re.compile(r"^contrast::([^:]+)::([^:]+)::([^:]+)$")
_SIMPLE_SPEC_RE = re.compile(
    r"^simple::([^:]+)::([^:]+)::([^:]+)\|within::([^:]+)::([^:]+)$"
)
_VALID_SIZE_FACTOR_TYPES = frozenset({"iterate", "poscounts", "ratio"})


def _prepare_count_table(table: biom.Table) -> pd.DataFrame:
    counts = table.to_dataframe(dense=True).round().astype(int)
    if (counts.values < 0).any():
        raise ValueError("Feature table counts must be non-negative integers.")
    return counts


def _filter_counts(counts: pd.DataFrame, min_total_count: int) -> pd.DataFrame:
    filtered = counts.loc[counts.sum(axis=1) >= min_total_count]
    if filtered.empty:
        raise ValueError(
            "No genes remain after filtering. Lower min_total_count or provide a denser feature table."
        )
    return filtered


def _collect_matching_samples(
    sample_ids: list[str], metadata_index: pd.Index, context: str
) -> list[str]:
    matched_samples = [sample_id for sample_id in sample_ids if sample_id in metadata_index]
    if len(matched_samples) < 2:
        raise ValueError(
            f"At least two samples must overlap between the feature table and {context}."
        )
    return matched_samples


def _normalize_formula(formula: str, parameter_name: str) -> str:
    normalized = str(formula).strip()
    if not normalized:
        raise ValueError(f"{parameter_name} is required.")
    if not _FORMULA_ALLOWED_CHARS.fullmatch(normalized):
        raise ValueError(
            f"{parameter_name} may only contain metadata column names, spaces, '+', ':', '*', and parentheses."
        )
    return normalized


def _extract_formula_columns(formula: str, parameter_name: str) -> list[str]:
    normalized = _normalize_formula(formula, parameter_name)
    columns = list(dict.fromkeys(_FORMULA_TOKEN_RE.findall(normalized)))
    if not columns:
        raise ValueError(
            f"{parameter_name} must reference at least one metadata column."
        )
    return columns


def _parse_reference_level_spec(reference_level_spec: str) -> tuple[str, str]:
    column, separator, level = str(reference_level_spec).strip().partition("::")
    if separator != "::" or not column or not level:
        raise ValueError(
            "reference_levels entries must use the form 'column::level'."
        )
    return column, level


def _normalize_reference_levels(reference_levels: list[str] | None) -> list[str]:
    normalized = []
    seen_columns = set()
    for raw_spec in reference_levels or []:
        spec = str(raw_spec).strip()
        if not spec:
            continue
        column, level = _parse_reference_level_spec(spec)
        if column in seen_columns:
            raise ValueError(
                f'Duplicate reference_levels entry provided for metadata column "{column}".'
            )
        normalized.append(f"{column}::{level}")
        seen_columns.add(column)
    return normalized


def _normalize_size_factor_type(size_factor_type: str) -> str:
    normalized = str(size_factor_type).strip().lower()
    if not normalized:
        raise ValueError("size_factor_type is required.")
    if normalized not in _VALID_SIZE_FACTOR_TYPES:
        supported_values = ", ".join(
            f'"{value}"' for value in sorted(_VALID_SIZE_FACTOR_TYPES)
        )
        raise ValueError(f"size_factor_type must be one of: {supported_values}.")
    return normalized


def _validate_effect_specs(effect_specs: list[str] | None, test: str) -> list[str]:
    normalized = []
    for raw_spec in effect_specs or []:
        spec = str(raw_spec).strip()
        if not spec:
            continue
        if (
            _COEF_SPEC_RE.match(spec) is None
            and _CONTRAST_SPEC_RE.match(spec) is None
            and _SIMPLE_SPEC_RE.match(spec) is None
        ):
            raise ValueError(
                "effect_specs entries must use one of: "
                "'coef::<resultsName>', "
                "'contrast::<factor>::<numerator>::<denominator>', or "
                "'simple::<factor>::<numerator>::<denominator>|within::<factor>::<level>'."
            )
        normalized.append(spec)

    if test == "lrt" and normalized:
        raise ValueError("effect_specs are not supported when test='lrt'.")

    return normalized


def _coerce_metadata_column(column: pd.Series) -> pd.Series:
    if is_numeric_dtype(column):
        return pd.to_numeric(column, errors="raise")
    return column.astype(str)


def _prepare_inputs(
    table: biom.Table,
    condition,
    min_total_count: int,
    reference_level: str,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], str]:
    counts = _prepare_count_table(table)
    sample_ids = list(table.ids(axis="sample"))

    metadata = condition.to_series().dropna().astype(str)
    matched_samples = _collect_matching_samples(
        sample_ids, metadata.index, "condition metadata"
    )

    counts = counts.loc[:, matched_samples]
    metadata = metadata.loc[matched_samples]

    if metadata.nunique() < 2:
        raise ValueError("Condition metadata must contain at least two unique levels.")

    counts = _filter_counts(counts, min_total_count)

    levels = sorted(metadata.unique().tolist())
    if not reference_level:
        if len(levels) != 2:
            raise ValueError(
                "Condition metadata has more than two levels. "
                "Set reference_level to define the baseline for all pairwise contrasts."
            )
        reference_level = levels[0]
    elif reference_level not in levels:
        raise ValueError(
            f'reference_level "{reference_level}" is not present in the condition metadata.'
        )

    comparison_levels = [level for level in levels if level != reference_level]
    if not comparison_levels:
        raise ValueError(
            "Condition metadata must contain at least one non-reference level."
        )

    coldata = pd.DataFrame({"condition": metadata})
    return counts, coldata, comparison_levels, reference_level


def _prepare_model_inputs(
    table: biom.Table,
    metadata,
    fixed_effects_formula: str,
    min_total_count: int,
    reference_levels: list[str] | None = None,
    reduced_formula: str = "",
) -> tuple[pd.DataFrame, pd.DataFrame, str, list[str], str]:
    normalized_formula = _normalize_formula(
        fixed_effects_formula, parameter_name="fixed_effects_formula"
    )
    normalized_reduced_formula = (
        _normalize_formula(reduced_formula, parameter_name="reduced_formula")
        if str(reduced_formula).strip()
        else ""
    )

    metadata_df = metadata.to_dataframe()
    sample_ids = list(table.ids(axis="sample"))
    counts = _prepare_count_table(table)
    referenced_columns = _extract_formula_columns(
        normalized_formula, parameter_name="fixed_effects_formula"
    )
    if normalized_reduced_formula:
        referenced_columns.extend(
            _extract_formula_columns(
                normalized_reduced_formula, parameter_name="reduced_formula"
            )
        )
    referenced_columns = list(dict.fromkeys(referenced_columns))

    missing_columns = [
        column for column in referenced_columns if column not in metadata_df.columns
    ]
    if missing_columns:
        missing_display = ", ".join(sorted(missing_columns))
        raise ValueError(
            "The metadata is missing columns required by the model formula: "
            f"{missing_display}."
        )

    normalized_reference_levels = _normalize_reference_levels(reference_levels)
    for reference_level_spec in normalized_reference_levels:
        column, _ = _parse_reference_level_spec(reference_level_spec)
        if column not in referenced_columns:
            raise ValueError(
                f'reference_levels entry "{reference_level_spec}" refers to "{column}", '
                "but that column is not used in the fitted model."
            )

    matched_samples = _collect_matching_samples(
        sample_ids, metadata_df.index, "metadata"
    )
    coldata = metadata_df.loc[matched_samples, referenced_columns].copy()
    counts = counts.loc[:, matched_samples]

    coldata = coldata.loc[~coldata.isna().any(axis=1)].copy()
    if coldata.shape[0] < 2:
        raise ValueError(
            "At least two samples with complete metadata are required after dropping missing values."
        )

    counts = counts.loc[:, coldata.index]
    counts = _filter_counts(counts, min_total_count)
    for column in coldata.columns:
        coldata[column] = _coerce_metadata_column(coldata[column])

    for reference_level_spec in normalized_reference_levels:
        column, level = _parse_reference_level_spec(reference_level_spec)
        if is_numeric_dtype(coldata[column]):
            raise ValueError(
                f'reference_levels entry "{reference_level_spec}" refers to numeric metadata column "{column}".'
            )
        observed_levels = set(coldata[column].astype(str))
        if level not in observed_levels:
            raise ValueError(
                f'reference_levels entry "{reference_level_spec}" refers to level "{level}", '
                f'which is not present in metadata column "{column}".'
            )

    return (
        counts,
        coldata,
        normalized_formula,
        normalized_reference_levels,
        normalized_reduced_formula,
    )


def _compute_normalized_counts(
    counts_df: pd.DataFrame, size_factors: pd.Series
) -> pd.DataFrame:
    aligned_size_factors = size_factors.reindex(counts_df.columns.map(str))
    missing_samples = aligned_size_factors[aligned_size_factors.isna()].index.tolist()
    if missing_samples:
        missing_display = ", ".join(missing_samples)
        raise RuntimeError(
            "DESeq2 completed but did not return size factors for sample(s): "
            f"{missing_display}."
        )

    normalized = counts_df.div(aligned_size_factors.astype(float), axis=1)
    normalized.index = normalized.index.map(str)
    normalized.columns = normalized.columns.map(str)
    normalized.insert(0, "feature_id", normalized.index)
    normalized = normalized.reset_index(drop=True)
    normalized.columns.name = None
    return normalized


def _compute_sample_distance_matrix(vst_counts: pd.DataFrame) -> pd.DataFrame:
    sample_ids = vst_counts.columns.map(str)
    sample_matrix = vst_counts.to_numpy(dtype=float).T
    squared_norms = np.sum(sample_matrix * sample_matrix, axis=1, keepdims=True)
    squared_distances = squared_norms + squared_norms.T - 2 * (
        sample_matrix @ sample_matrix.T
    )
    squared_distances = np.maximum(squared_distances, 0.0)
    distances = np.sqrt(squared_distances)
    return pd.DataFrame(distances, index=sample_ids, columns=sample_ids)


def _compute_sample_distance_order(
    sample_distance_matrix: pd.DataFrame,
) -> tuple[str, ...]:
    sample_ids = tuple(sample_distance_matrix.index.map(str))
    sample_count = len(sample_ids)
    if sample_count < 2:
        return sample_ids

    distances = sample_distance_matrix.to_numpy(dtype=float)
    clusters = {
        index: {"order": (index,), "height": 0.0}
        for index in range(sample_count)
    }
    cluster_distances = {
        (left, right): float(distances[left, right])
        for left in range(sample_count)
        for right in range(left + 1, sample_count)
    }

    active_clusters = list(range(sample_count))
    next_cluster_id = sample_count

    while len(active_clusters) > 1:
        best_pair = None
        best_key = None
        active_clusters_sorted = sorted(
            active_clusters, key=lambda cluster_id: clusters[cluster_id]["order"]
        )
        for left_index, left_cluster in enumerate(active_clusters_sorted[:-1]):
            for right_cluster in active_clusters_sorted[left_index + 1 :]:
                pair = (
                    (left_cluster, right_cluster)
                    if left_cluster < right_cluster
                    else (right_cluster, left_cluster)
                )
                candidate_key = (
                    cluster_distances[pair],
                    clusters[left_cluster]["order"],
                    clusters[right_cluster]["order"],
                )
                if best_key is None or candidate_key < best_key:
                    best_key = candidate_key
                    best_pair = pair

        if best_pair is None or best_key is None:
            break

        left_cluster, right_cluster = sorted(
            best_pair,
            key=lambda cluster_id: (
                clusters[cluster_id]["height"],
                clusters[cluster_id]["order"],
            ),
        )
        merged_order = (
            clusters[left_cluster]["order"] + clusters[right_cluster]["order"]
        )
        merged_height = float(best_key[0])
        clusters[next_cluster_id] = {
            "order": merged_order,
            "height": merged_height,
        }

        remaining_clusters = [
            cluster_id
            for cluster_id in active_clusters
            if cluster_id not in best_pair
        ]
        for other_cluster in remaining_clusters:
            left_distance = cluster_distances[
                (
                    (left_cluster, other_cluster)
                    if left_cluster < other_cluster
                    else (other_cluster, left_cluster)
                )
            ]
            right_distance = cluster_distances[
                (
                    (right_cluster, other_cluster)
                    if right_cluster < other_cluster
                    else (other_cluster, right_cluster)
                )
            ]
            pair = (
                (other_cluster, next_cluster_id)
                if other_cluster < next_cluster_id
                else (next_cluster_id, other_cluster)
            )
            cluster_distances[pair] = max(left_distance, right_distance)

        active_clusters = remaining_clusters + [next_cluster_id]
        next_cluster_id += 1

    final_order = clusters[active_clusters[0]]["order"]
    return tuple(sample_ids[index] for index in final_order)


def _compute_sample_pca(
    vst_counts: pd.DataFrame,
) -> tuple[pd.DataFrame, tuple[float, float]]:
    sample_ids = vst_counts.columns.map(str)
    sample_count = len(sample_ids)
    if sample_count == 0:
        return pd.DataFrame(columns=["PC1", "PC2"]), ()

    ntop = min(500, vst_counts.shape[0])
    if ntop == 0:
        empty_scores = pd.DataFrame(
            {"PC1": np.zeros(sample_count), "PC2": np.zeros(sample_count)},
            index=sample_ids,
        )
        empty_scores.index.name = None
        return empty_scores, ()

    matrix = vst_counts.to_numpy(dtype=float)
    ddof = 1 if sample_count > 1 else 0
    row_variances = np.var(matrix, axis=1, ddof=ddof)
    selected_rows = np.argsort(-row_variances, kind="mergesort")[:ntop]
    selected_matrix = matrix[selected_rows, :].T
    centered_matrix = selected_matrix - selected_matrix.mean(axis=0, keepdims=True)

    if centered_matrix.size == 0:
        singular_values = np.array([], dtype=float)
        sample_scores = np.zeros((sample_count, 0), dtype=float)
    else:
        left_singular_vectors, singular_values, _ = np.linalg.svd(
            centered_matrix, full_matrices=False
        )
        sample_scores = left_singular_vectors * singular_values

    explained_variance = singular_values * singular_values
    total_variance = float(explained_variance.sum())
    if total_variance > 0:
        percent_variance = explained_variance / total_variance
    else:
        percent_variance = np.zeros(max(1, len(singular_values)), dtype=float)
        percent_variance[0] = 1.0

    pc1_values = (
        sample_scores[:, 0]
        if sample_scores.shape[1] >= 1
        else np.zeros(sample_count, dtype=float)
    )
    pc2_values = (
        sample_scores[:, 1]
        if sample_scores.shape[1] >= 2
        else np.zeros(sample_count, dtype=float)
    )
    percent_pc1 = (
        float(percent_variance[0] * 100) if len(percent_variance) >= 1 else 0.0
    )
    percent_pc2 = (
        float(percent_variance[1] * 100) if len(percent_variance) >= 2 else 0.0
    )

    sample_pca_scores = pd.DataFrame(
        {"PC1": pc1_values, "PC2": pc2_values},
        index=sample_ids,
    )
    sample_pca_scores.index.name = None
    return sample_pca_scores, (percent_pc1, percent_pc2)


def _compute_count_matrix_heatmap(
    vst_counts: pd.DataFrame,
    normalized_counts: pd.DataFrame,
    sample_distance_order: tuple[str, ...],
) -> pd.DataFrame:
    ordered_sample_ids = [
        sample_id
        for sample_id in sample_distance_order
        if sample_id in vst_counts.columns
    ]
    if not ordered_sample_ids:
        ordered_sample_ids = list(vst_counts.columns.map(str))

    normalized_by_feature = normalized_counts.set_index("feature_id")
    feature_means = normalized_by_feature.mean(axis=1)
    top_heatmap_n = min(100, len(feature_means))
    selected_feature_ids = feature_means.sort_values(
        ascending=False, kind="mergesort"
    ).index[:top_heatmap_n]
    heatmap = vst_counts.loc[selected_feature_ids, ordered_sample_ids].copy()
    heatmap.index = heatmap.index.map(str)
    heatmap.columns = heatmap.columns.map(str)
    heatmap.index.name = None
    heatmap.columns.name = None
    return heatmap


def _compute_run_analytics(
    counts_df: pd.DataFrame,
    size_factors: pd.Series,
    vst_counts: pd.DataFrame,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    tuple[str, ...],
    pd.DataFrame,
    tuple[float, float],
    pd.DataFrame,
]:
    expected_sample_ids = list(counts_df.columns.map(str))
    missing_vst_samples = [
        sample_id for sample_id in expected_sample_ids if sample_id not in vst_counts.columns
    ]
    if missing_vst_samples:
        missing_display = ", ".join(missing_vst_samples)
        raise RuntimeError(
            "DESeq2 completed but the VST matrix is missing sample(s): "
            f"{missing_display}."
        )

    expected_feature_ids = list(counts_df.index.map(str))
    missing_vst_features = [
        feature_id for feature_id in expected_feature_ids if feature_id not in vst_counts.index
    ]
    if missing_vst_features:
        missing_display = ", ".join(missing_vst_features[:10])
        if len(missing_vst_features) > 10:
            missing_display += ", ..."
        raise RuntimeError(
            "DESeq2 completed but the VST matrix is missing feature(s): "
            f"{missing_display}."
        )

    vst_counts = vst_counts.loc[expected_feature_ids, expected_sample_ids]
    normalized_counts = _compute_normalized_counts(counts_df, size_factors)
    sample_distance_matrix = _compute_sample_distance_matrix(vst_counts)
    sample_distance_matrix.index.name = None
    sample_distance_matrix.columns.name = None
    sample_distance_order = _compute_sample_distance_order(sample_distance_matrix)
    sample_pca_scores, sample_pca_percent_variance = _compute_sample_pca(vst_counts)
    count_matrix_heatmap = _compute_count_matrix_heatmap(
        vst_counts,
        normalized_counts,
        sample_distance_order,
    )
    return (
        normalized_counts,
        sample_distance_matrix,
        sample_distance_order,
        sample_pca_scores,
        sample_pca_percent_variance,
        count_matrix_heatmap,
    )


def _write_lines(path: Path, lines: list[str]) -> None:
    payload = "\n".join(lines)
    if payload:
        payload += "\n"
    path.write_text(payload, encoding="utf-8")


def _first_non_empty_string(value) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return text


def _first_value_from_column(frame: pd.DataFrame, column_name: str) -> str:
    if column_name not in frame.columns:
        return ""
    for value in frame[column_name]:
        text = _first_non_empty_string(value)
        if text:
            return text
    return ""


def _unique_non_empty_values(values) -> tuple[str, ...]:
    normalized = []
    seen = set()
    for value in values:
        text = _first_non_empty_string(value)
        if text and text not in seen:
            normalized.append(text)
            seen.add(text)
    return tuple(normalized)


def _run_deseq2_with_frames(
    counts_df: pd.DataFrame,
    coldata_df: pd.DataFrame,
    fixed_effects_formula: str,
    reference_levels: list[str] | None = None,
    effect_specs: list[str] | None = None,
    test: str = "wald",
    reduced_formula: str = "",
    fit_type: str = "parametric",
    size_factor_type: str = "poscounts",
    alpha: float = 0.05,
    cooks_cutoff: bool = True,
    independent_filtering: bool = True,
    legacy_test_level: str = "",
    legacy_reference_level: str = "",
) -> DESeq2RunResult:
    normalized_reference_levels = _normalize_reference_levels(reference_levels)
    normalized_size_factor_type = _normalize_size_factor_type(size_factor_type)
    normalized_test = str(test).strip().lower() or "wald"
    if normalized_test not in {"wald", "lrt"}:
        raise ValueError("test must be either 'wald' or 'lrt'.")
    normalized_effect_specs = _validate_effect_specs(effect_specs, normalized_test)
    normalized_reduced_formula = (
        _normalize_formula(reduced_formula, parameter_name="reduced_formula")
        if str(reduced_formula).strip()
        else ""
    )
    if normalized_test == "lrt" and not normalized_reduced_formula:
        raise ValueError("reduced_formula is required when test='lrt'.")

    with tempfile.TemporaryDirectory(prefix="q2-deseq2-") as temp_dir:
        temp_path = Path(temp_dir)
        counts_fp = temp_path / "counts.tsv"
        coldata_fp = temp_path / "coldata.tsv"
        results_fp = temp_path / "deseq2_results.tsv"
        size_factors_fp = temp_path / "size_factors.tsv"
        vst_counts_fp = temp_path / "vst_counts.tsv"
        summary_fp = temp_path / "deseq2_summary.txt"
        results_names_fp = temp_path / "results_names.txt"
        reference_levels_fp = temp_path / "reference_levels.txt"
        effect_specs_fp = temp_path / "effect_specs.txt"
        script_fp = Path(__file__).parent / "run_deseq2.R"
        ma_plot_fp = temp_path / "ma_plot.png"

        counts_df.to_csv(counts_fp, sep="\t", index_label="feature_id")
        coldata_df.to_csv(coldata_fp, sep="\t", index_label="sample_id")
        _write_lines(reference_levels_fp, normalized_reference_levels)
        _write_lines(effect_specs_fp, normalized_effect_specs)

        cmd = [
            "Rscript",
            str(script_fp),
            "--counts",
            str(counts_fp),
            "--coldata",
            str(coldata_fp),
            "--results",
            str(results_fp),
            "--size-factors",
            str(size_factors_fp),
            "--vst-counts",
            str(vst_counts_fp),
            "--summary",
            str(summary_fp),
            "--ma-plot",
            str(ma_plot_fp),
            "--results-names",
            str(results_names_fp),
            "--reference-levels",
            str(reference_levels_fp),
            "--effect-specs",
            str(effect_specs_fp),
            "--fit-type",
            fit_type,
            "--size-factor-type",
            normalized_size_factor_type,
            "--alpha",
            str(alpha),
            "--cooks-cutoff",
            str(cooks_cutoff).lower(),
            "--independent-filtering",
            str(independent_filtering).lower(),
            "--fixed-effects-formula",
            fixed_effects_formula,
            "--test",
            normalized_test,
            "--reduced-formula",
            normalized_reduced_formula,
        ]

        try:
            run(cmd, check=True, capture_output=True, text=True)
        except CalledProcessError as exc:
            detail = exc.stderr.strip() or exc.stdout.strip()
            raise RuntimeError(
                f"DESeq2 command failed with exit code {exc.returncode}: {detail}"
            ) from exc

        expected_outputs = {
            "deseq2_results.tsv": results_fp,
            "size_factors.tsv": size_factors_fp,
            "vst_counts.tsv": vst_counts_fp,
            "results_names.txt": results_names_fp,
        }
        for expected_name, path in expected_outputs.items():
            if not path.exists():
                raise RuntimeError(
                    f'DESeq2 completed but expected output file "{expected_name}" was not created.'
                )

        results_df = pd.read_csv(results_fp, sep="\t")
        size_factors_df = pd.read_csv(size_factors_fp, sep="\t")
        required_size_factor_columns = {"sample_id", "size_factor"}
        if not required_size_factor_columns <= set(size_factors_df.columns):
            missing_columns = ", ".join(
                sorted(required_size_factor_columns - set(size_factors_df.columns))
            )
            raise RuntimeError(
                "DESeq2 completed but size_factors.tsv is missing required column(s): "
                f"{missing_columns}."
            )
        size_factors = pd.Series(
            size_factors_df["size_factor"].astype(float).to_numpy(),
            index=size_factors_df["sample_id"].astype(str),
        )
        vst_counts = pd.read_csv(vst_counts_fp, sep="\t", index_col=0)
        vst_counts.index = vst_counts.index.map(str)
        vst_counts.columns = vst_counts.columns.map(str)
        vst_counts.index.name = None
        vst_counts.columns.name = None
        (
            normalized_counts_df,
            sample_distance_matrix,
            sample_distance_order,
            sample_pca_scores,
            sample_pca_percent_variance,
            count_matrix_heatmap,
        ) = _compute_run_analytics(
            counts_df=counts_df,
            size_factors=size_factors,
            vst_counts=vst_counts,
        )
        available_results_names = _unique_non_empty_values(
            results_names_fp.read_text(encoding="utf-8").splitlines()
        )
        selected_effect_specs = _unique_non_empty_values(
            results_df["effect_id"] if "effect_id" in results_df.columns else ()
        )
        default_effect_id = _first_value_from_column(results_df, "effect_id")
        sample_metadata_df = coldata_df.copy()
        sample_metadata_df.index = sample_metadata_df.index.map(str)
        sample_metadata_df.columns = sample_metadata_df.columns.map(str)

        return DESeq2RunResult(
            results=results_df,
            normalized_counts=normalized_counts_df,
            sample_metadata=sample_metadata_df,
            test_level=legacy_test_level
            or _first_value_from_column(results_df, "test_level"),
            reference_level=legacy_reference_level
            or _first_value_from_column(results_df, "reference_level"),
            default_effect_id=default_effect_id,
            fixed_effects_formula=fixed_effects_formula,
            reference_levels=tuple(normalized_reference_levels),
            test=normalized_test,
            reduced_formula=normalized_reduced_formula,
            available_results_names=available_results_names,
            selected_effect_specs=selected_effect_specs,
            sample_distance_matrix=sample_distance_matrix,
            sample_distance_order=sample_distance_order,
            sample_pca_scores=sample_pca_scores,
            sample_pca_percent_variance=sample_pca_percent_variance,
            count_matrix_heatmap=count_matrix_heatmap,
        )


def run_deseq2_model(
    table: biom.Table,
    metadata,
    fixed_effects_formula: str,
    reference_levels: list[str] | None = None,
    effect_specs: list[str] | None = None,
    test: str = "wald",
    reduced_formula: str = "",
    min_total_count: int = 10,
    fit_type: str = "parametric",
    size_factor_type: str = "poscounts",
    alpha: float = 0.05,
    cooks_cutoff: bool = True,
    independent_filtering: bool = True,
) -> DESeq2RunResult:
    (
        counts_df,
        coldata_df,
        normalized_formula,
        normalized_reference_levels,
        normalized_reduced_formula,
    ) = _prepare_model_inputs(
        table=table,
        metadata=metadata,
        fixed_effects_formula=fixed_effects_formula,
        min_total_count=min_total_count,
        reference_levels=reference_levels,
        reduced_formula=reduced_formula,
    )

    return _run_deseq2_with_frames(
        counts_df=counts_df,
        coldata_df=coldata_df,
        fixed_effects_formula=normalized_formula,
        reference_levels=normalized_reference_levels,
        effect_specs=effect_specs,
        test=test,
        reduced_formula=normalized_reduced_formula,
        fit_type=fit_type,
        size_factor_type=size_factor_type,
        alpha=alpha,
        cooks_cutoff=cooks_cutoff,
        independent_filtering=independent_filtering,
    )


def run_deseq2(
    table: biom.Table,
    condition,
    reference_level: str = "",
    min_total_count: int = 10,
    fit_type: str = "parametric",
    size_factor_type: str = "poscounts",
    alpha: float = 0.05,
    cooks_cutoff: bool = True,
    independent_filtering: bool = True,
) -> DESeq2RunResult:
    counts_df, coldata_df, comparison_levels, reference_level = _prepare_inputs(
        table, condition, min_total_count, reference_level
    )
    effect_specs = [
        f"contrast::condition::{comparison_level}::{reference_level}"
        for comparison_level in comparison_levels
    ]

    return _run_deseq2_with_frames(
        counts_df=counts_df,
        coldata_df=coldata_df,
        fixed_effects_formula="condition",
        reference_levels=[f"condition::{reference_level}"],
        effect_specs=effect_specs,
        test="wald",
        fit_type=fit_type,
        size_factor_type=size_factor_type,
        alpha=alpha,
        cooks_cutoff=cooks_cutoff,
        independent_filtering=independent_filtering,
        legacy_test_level=comparison_levels[0],
        legacy_reference_level=reference_level,
    )


def _estimate_model(
    table: biom.Table,
    metadata: str,
    fixed_effects_formula: str,
    reference_levels: list[str] | None = None,
    effect_specs: list[str] | None = None,
    test: str = "wald",
    reduced_formula: str = "",
    min_total_count: int = 10,
    fit_type: str = "parametric",
    size_factor_type: str = "poscounts",
    alpha: float = 0.05,
    cooks_cutoff: bool = True,
    independent_filtering: bool = True,
) -> (pd.DataFrame, DESeq2RunDirectoryFormat):
    run_result = run_deseq2_model(
        table=table,
        metadata=metadata,
        fixed_effects_formula=fixed_effects_formula,
        reference_levels=reference_levels,
        effect_specs=effect_specs,
        test=test,
        reduced_formula=reduced_formula,
        min_total_count=min_total_count,
        fit_type=fit_type,
        size_factor_type=size_factor_type,
        alpha=alpha,
        cooks_cutoff=cooks_cutoff,
        independent_filtering=independent_filtering,
    )
    run_data = write_run_result_artifact(run_result=run_result, alpha=alpha)
    return run_result.results, run_data


def _estimate_differential_expression(
    table: biom.Table,
    condition: str,
    reference_level: str = "",
    min_total_count: int = 10,
    fit_type: str = "parametric",
    size_factor_type: str = "poscounts",
    alpha: float = 0.05,
    cooks_cutoff: bool = True,
    independent_filtering: bool = True,
) -> (pd.DataFrame, DESeq2RunDirectoryFormat):
    run_result = run_deseq2(
        table=table,
        condition=condition,
        reference_level=reference_level,
        min_total_count=min_total_count,
        fit_type=fit_type,
        size_factor_type=size_factor_type,
        alpha=alpha,
        cooks_cutoff=cooks_cutoff,
        independent_filtering=independent_filtering,
    )
    run_data = write_run_result_artifact(run_result=run_result, alpha=alpha)
    return run_result.results, run_data
