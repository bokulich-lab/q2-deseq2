# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Lesser General Public License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
from pathlib import Path
from shutil import copytree
from urllib.parse import unquote

import pandas as pd
import q2templates
from pandas.api.types import is_numeric_dtype
from q2_types.genome_data import LociDirectoryFormat

from q2_deseq2.types import DESeq2RunDirectoryFormat
from q2_deseq2.utils.run_data import DESeq2RunResult, _parse_run_results

_ASSETS_DIR = Path(__file__).resolve().parent / "assets"


def _value_or_none(value):
    if pd.isna(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _parse_gff3_attributes(raw_attributes: str) -> dict[str, str]:
    parsed = {}
    for entry in raw_attributes.split(";"):
        entry = entry.strip()
        if not entry:
            continue
        if "=" in entry:
            key, value = entry.split("=", 1)
        else:
            key, value = entry, ""
        parsed[key.strip()] = unquote(value.strip())
    return parsed


def _first_non_empty(attributes: dict[str, str], keys: list[str]) -> str | None:
    for key in keys:
        value = attributes.get(key, "").strip()
        if value:
            return value.split(",")[0]
    return None


def _resolve_reference_gff(annotations: LociDirectoryFormat, reference_id: str) -> Path:
    if not reference_id:
        raise ValueError("reference_id is required when gene annotations are provided.")

    file_map = annotations.file_dict()
    gff_fp = file_map.get(reference_id)
    if gff_fp is None:
        available_refs = ", ".join(sorted(file_map))
        raise ValueError(
            f'Reference "{reference_id}" was not found in annotations artifact. '
            f"Available references: {available_refs}"
        )

    return Path(gff_fp)


def _load_loci(gff_path: Path) -> pd.DataFrame:
    annotation_lookup = {}
    with gff_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            attributes = _parse_gff3_attributes(parts[8])
            feature_ids = []
            for key in ["ID", "locus_tag", "Name", "gene", "Parent"]:
                value = attributes.get(key, "").strip()
                if value:
                    feature_ids.extend(v.strip() for v in value.split(",") if v.strip())

            if not feature_ids:
                continue

            gene_name = _first_non_empty(
                attributes, ["gene", "Name", "gene_name", "locus_tag"]
            )
            product = _first_non_empty(attributes, ["product", "description", "Note"])

            for feature_id in feature_ids:
                existing = annotation_lookup.get(
                    feature_id,
                    {"feature_id": feature_id, "gene_name": None, "product": None},
                )
                if existing["gene_name"] is None and gene_name is not None:
                    existing["gene_name"] = gene_name
                if existing["product"] is None and product is not None:
                    existing["product"] = product
                annotation_lookup[feature_id] = existing

    if not annotation_lookup:
        return pd.DataFrame(columns=["feature_id", "gene_name", "product"])

    annotation_table = pd.DataFrame(annotation_lookup.values())
    annotation_table["feature_id"] = annotation_table["feature_id"].astype(str)
    return annotation_table


def _add_annotations(
    result_table: pd.DataFrame, annotation_table: pd.DataFrame
) -> pd.DataFrame:
    merged = result_table.copy()
    merged["feature_id"] = merged["feature_id"].astype(str)
    if not annotation_table.empty:
        merged = merged.merge(annotation_table, how="left", on="feature_id")
    return merged


def _build_plot_records(result_table: pd.DataFrame) -> list[dict]:
    records = []
    for _, row in result_table.iterrows():
        feature_id = row.get("feature_id")
        if pd.isna(feature_id):
            feature_id = ""

        gene_name = row.get("gene_name")
        if pd.isna(gene_name):
            gene_name = None

        product = row.get("product")
        if pd.isna(product):
            product = None

        records.append(
            {
                "feature_id": str(feature_id),
                "gene_name": None if gene_name is None else str(gene_name),
                "product": None if product is None else str(product),
                "log2FoldChange": _value_or_none(row.get("log2FoldChange")),
                "pvalue": _value_or_none(row.get("pvalue")),
                "padj": _value_or_none(row.get("padj")),
                "baseMean": _value_or_none(row.get("baseMean")),
            }
        )

    return records


def _load_vega_spec(spec_name: str) -> dict:
    spec_path = _ASSETS_DIR / "vega" / f"{spec_name}.json"
    return json.loads(spec_path.read_text(encoding="utf-8"))


def _format_comparison_label(test_level: str, reference_level: str) -> str:
    return f"{test_level} vs. {reference_level}"


def _ensure_comparison_columns(
    result_table: pd.DataFrame,
    default_test_level: str,
    reference_level: str,
) -> pd.DataFrame:
    enriched = result_table.copy()
    if "test_level" not in enriched.columns:
        enriched["test_level"] = default_test_level
    else:
        enriched["test_level"] = enriched["test_level"].fillna(default_test_level)

    if "reference_level" not in enriched.columns:
        enriched["reference_level"] = reference_level
    else:
        enriched["reference_level"] = enriched["reference_level"].fillna(
            reference_level
        )

    if "comparison" not in enriched.columns:
        enriched["comparison"] = ""
    comparison_missing = enriched["comparison"].isna() | (
        enriched["comparison"].astype(str).str.strip() == ""
    )
    if comparison_missing.any():
        generated = enriched.loc[comparison_missing].apply(
            lambda row: _format_comparison_label(
                str(row["test_level"]), str(row["reference_level"])
            ),
            axis=1,
        )
        enriched["comparison"] = enriched["comparison"].astype(object)
        enriched.loc[comparison_missing, "comparison"] = generated.values

    return enriched


def _legacy_effect_id(comparison: str, test_level: str, reference_level: str) -> str:
    if test_level and reference_level:
        return f"legacy::{test_level}::{reference_level}"
    if comparison:
        return f"legacy::{comparison}"
    return "legacy::effect"


def _ensure_effect_columns(
    result_table: pd.DataFrame,
    default_effect_id: str,
    default_test_level: str,
    reference_level: str,
) -> tuple[pd.DataFrame, str]:
    enriched = _ensure_comparison_columns(
        result_table,
        default_test_level=default_test_level,
        reference_level=reference_level,
    )

    if "effect_id" not in enriched.columns:
        enriched["effect_id"] = ""
    if "effect_label" not in enriched.columns:
        enriched["effect_label"] = ""
    if "effect_kind" not in enriched.columns:
        enriched["effect_kind"] = ""
    if "effect_expression" not in enriched.columns:
        enriched["effect_expression"] = ""

    for index, row in enriched.iterrows():
        comparison = str(row.get("comparison", "") or "").strip()
        test_level = str(row.get("test_level", "") or "").strip()
        reference = str(row.get("reference_level", "") or "").strip()
        effect_id = str(row.get("effect_id", "") or "").strip()
        effect_label = str(row.get("effect_label", "") or "").strip()
        effect_kind = str(row.get("effect_kind", "") or "").strip()
        effect_expression = str(row.get("effect_expression", "") or "").strip()

        if not effect_id:
            effect_id = _legacy_effect_id(comparison, test_level, reference)
            enriched.at[index, "effect_id"] = effect_id
        if not effect_label:
            effect_label = comparison or effect_id
            enriched.at[index, "effect_label"] = effect_label
        if not effect_kind:
            enriched.at[index, "effect_kind"] = "comparison"
        if not effect_expression:
            enriched.at[index, "effect_expression"] = comparison or effect_label

    available_effect_ids = [
        str(effect_id).strip()
        for effect_id in enriched["effect_id"].dropna().tolist()
        if str(effect_id).strip()
    ]
    if default_effect_id not in available_effect_ids:
        default_effect_id = available_effect_ids[0] if available_effect_ids else ""

    return enriched, default_effect_id


def _collect_effect_options(
    result_table: pd.DataFrame, default_effect_id: str
) -> tuple[list[dict[str, str]], str, str]:
    effect_frame = (
        result_table.loc[:, ["effect_id", "effect_label", "effect_kind"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    effect_options = effect_frame.to_dict(orient="records")
    if not effect_options:
        effect_options = [
            {
                "effect_id": default_effect_id or "legacy::effect",
                "effect_label": "Effect",
                "effect_kind": "comparison",
            }
        ]

    known_effect_ids = {option["effect_id"] for option in effect_options}
    if default_effect_id not in known_effect_ids:
        default_effect_id = effect_options[0]["effect_id"]

    effect_options.sort(key=lambda option: option["effect_id"] != default_effect_id)
    default_effect_label = next(
        option["effect_label"]
        for option in effect_options
        if option["effect_id"] == default_effect_id
    )
    return effect_options, default_effect_id, default_effect_label


def _summarize_results(
    result_table: pd.DataFrame, alpha: float
) -> dict[str, int | str]:
    log2fc = pd.to_numeric(result_table.get("log2FoldChange"), errors="coerce")
    padj = pd.to_numeric(result_table.get("padj"), errors="coerce")
    plottable_mask = log2fc.notna() & padj.notna() & (padj > 0)
    significant_mask = plottable_mask & (padj <= alpha)

    summary = {
        "total_features": int(len(result_table)),
        "plottable_features": int(plottable_mask.sum()),
        "significant_features": int(significant_mask.sum()),
        "alpha_display": str(alpha),
    }
    if {"gene_name", "product"} & set(result_table.columns):
        annotated_mask = pd.Series(False, index=result_table.index)
        for column in ["gene_name", "product"]:
            if column in result_table.columns:
                annotated_mask = annotated_mask | result_table[column].notna()
        summary["annotated_features"] = int(annotated_mask.sum())

    return summary


def _prepare_table_payload(result_table: pd.DataFrame) -> str:
    serializable = result_table.astype(object).where(pd.notna(result_table), None)
    columns = [str(column) for column in serializable.columns]
    rows = serializable.loc[:, columns].values.tolist()
    return json.dumps({"columns": columns, "data": rows}).replace("NaN", "null")


def _resolve_sample_distance_order(
    sample_distance_matrix: pd.DataFrame, sample_distance_order: tuple[str, ...]
) -> list[str]:
    available_ids = [
        str(sample_id) for sample_id in sample_distance_matrix.index.tolist()
    ]
    ordered_ids = []
    seen = set()
    for sample_id in sample_distance_order:
        normalized = str(sample_id).strip()
        if normalized and normalized in available_ids and normalized not in seen:
            ordered_ids.append(normalized)
            seen.add(normalized)

    for sample_id in available_ids:
        if sample_id not in seen:
            ordered_ids.append(sample_id)

    return ordered_ids


def _prepare_sample_distance_payload(
    sample_distance_matrix: pd.DataFrame,
    sample_distance_order: tuple[str, ...],
    sample_metadata: pd.DataFrame | None = None,
    reference_levels: tuple[str, ...] = (),
) -> str:
    matrix = sample_distance_matrix.copy()
    matrix.index = matrix.index.map(str)
    matrix.columns = matrix.columns.map(str)
    ordered_ids = _resolve_sample_distance_order(matrix, sample_distance_order)
    sample_labels, sample_metadata_text = _build_sample_label_map(
        sample_metadata=sample_metadata,
        ordered_sample_ids=ordered_ids,
        reference_levels=reference_levels,
    )

    records = []
    for sample_y in ordered_ids:
        if sample_y not in matrix.index:
            continue
        for sample_x in ordered_ids:
            if sample_x not in matrix.columns:
                continue
            records.append(
                {
                    "sample_x": sample_x,
                    "sample_y": sample_y,
                    "sample_y_label": sample_labels.get(sample_y, sample_y),
                    "distance": _value_or_none(matrix.at[sample_y, sample_x]),
                    "is_diagonal": sample_x == sample_y,
                    "sample_x_metadata": sample_metadata_text.get(sample_x, ""),
                    "sample_y_metadata": sample_metadata_text.get(sample_y, ""),
                }
            )

    return json.dumps(records).replace("NaN", "null")


def _reference_level_map(reference_levels: tuple[str, ...]) -> dict[str, str]:
    mapping = {}
    for raw_spec in reference_levels:
        column, separator, level = str(raw_spec).strip().partition("::")
        if separator == "::" and column and level and column not in mapping:
            mapping[column] = level
    return mapping


def _prepare_sample_metadata_frame(
    sample_metadata: pd.DataFrame | None, ordered_sample_ids: list[str]
) -> pd.DataFrame | None:
    if sample_metadata is None or sample_metadata.empty:
        return None

    metadata = sample_metadata.copy()
    metadata.index = metadata.index.map(str)
    metadata.columns = metadata.columns.map(str)

    available_sample_ids = [
        sample_id for sample_id in ordered_sample_ids if sample_id in metadata.index
    ]
    if not available_sample_ids:
        available_sample_ids = metadata.index.tolist()

    metadata = metadata.loc[available_sample_ids]
    return metadata


def _ordered_sample_metadata_columns(
    sample_metadata: pd.DataFrame | None, reference_levels: tuple[str, ...]
) -> list[str]:
    if sample_metadata is None or sample_metadata.empty:
        return []

    reference_map = _reference_level_map(reference_levels)
    candidate_columns = []
    for column in sample_metadata.columns:
        series = sample_metadata[column]
        non_missing = series.dropna()
        if non_missing.empty:
            continue
        if is_numeric_dtype(series) and column not in reference_map:
            continue
        candidate_columns.append(column)

    ordered_columns = []
    seen_columns = set()
    for column in list(reference_map) + candidate_columns:
        if column in candidate_columns and column not in seen_columns:
            ordered_columns.append(column)
            seen_columns.add(column)

    return ordered_columns


def _format_metadata_column_label(
    column: str, reference_levels: tuple[str, ...]
) -> str:
    reference_level = _reference_level_map(reference_levels).get(column)
    if reference_level:
        return f"{column} (ref: {reference_level})"
    return column


def _build_sample_label_map(
    sample_metadata: pd.DataFrame | None,
    ordered_sample_ids: list[str],
    reference_levels: tuple[str, ...],
) -> tuple[dict[str, str], dict[str, str]]:
    metadata = _prepare_sample_metadata_frame(sample_metadata, ordered_sample_ids)
    if metadata is None or metadata.empty:
        return {}, {}

    available_sample_ids = metadata.index.tolist()
    ordered_columns = _ordered_sample_metadata_columns(metadata, reference_levels)

    sample_labels = {}
    sample_metadata_text = {}
    for sample_id in available_sample_ids:
        parts = []
        for column in ordered_columns:
            value = metadata.at[sample_id, column]
            if pd.isna(value):
                continue
            text = str(value).strip()
            if not text:
                continue
            parts.append(f"{column}={text}")

        sample_metadata_text[sample_id] = "; ".join(parts)
        sample_labels[sample_id] = (
            sample_id if not parts else f"{sample_id} | " + " | ".join(parts)
        )

    return sample_labels, sample_metadata_text


def _prepare_sample_pca_payload(
    sample_pca_scores: pd.DataFrame,
    sample_pca_percent_variance: tuple[float, float] = (),
    sample_distance_order: tuple[str, ...] = (),
    sample_metadata: pd.DataFrame | None = None,
    reference_levels: tuple[str, ...] = (),
) -> str:
    scores = sample_pca_scores.copy()
    scores.index = scores.index.map(str)
    scores.columns = scores.columns.map(str)

    ordered_ids = []
    seen = set()
    for sample_id in sample_distance_order:
        normalized = str(sample_id).strip()
        if normalized and normalized in scores.index and normalized not in seen:
            ordered_ids.append(normalized)
            seen.add(normalized)
    for sample_id in scores.index.tolist():
        if sample_id not in seen:
            ordered_ids.append(sample_id)

    metadata = _prepare_sample_metadata_frame(sample_metadata, ordered_ids)
    ordered_columns = _ordered_sample_metadata_columns(metadata, reference_levels)
    group_field = ordered_columns[0] if ordered_columns else ""
    group_label = (
        _format_metadata_column_label(group_field, reference_levels)
        if group_field
        else ""
    )
    sample_labels, sample_metadata_text = _build_sample_label_map(
        sample_metadata=metadata,
        ordered_sample_ids=ordered_ids,
        reference_levels=reference_levels,
    )

    points = []
    for sample_id in ordered_ids:
        if sample_id not in scores.index:
            continue

        pc1 = _value_or_none(scores.at[sample_id, "PC1"])
        pc2 = _value_or_none(scores.at[sample_id, "PC2"])
        if pc1 is None or pc2 is None:
            continue

        group_value = ""
        if metadata is not None and group_field and sample_id in metadata.index:
            raw_group_value = metadata.at[sample_id, group_field]
            if not pd.isna(raw_group_value):
                group_value = str(raw_group_value).strip()

        points.append(
            {
                "sample_id": sample_id,
                "sample_label": sample_labels.get(sample_id, sample_id),
                "sample_metadata": sample_metadata_text.get(sample_id, ""),
                "group_value": group_value or "Samples",
                "PC1": pc1,
                "PC2": pc2,
            }
        )

    percent_pc1 = (
        float(sample_pca_percent_variance[0])
        if len(sample_pca_percent_variance) >= 1
        else 0.0
    )
    percent_pc2 = (
        float(sample_pca_percent_variance[1])
        if len(sample_pca_percent_variance) >= 2
        else 0.0
    )

    payload = {
        "points": points,
        "group_field": group_field,
        "group_label": group_label,
        "percent_variance": {
            "PC1": percent_pc1,
            "PC2": percent_pc2,
        },
    }
    return json.dumps(payload).replace("NaN", "null")


def _prepare_count_matrix_heatmap_payload(
    count_matrix_heatmap: pd.DataFrame,
    sample_metadata: pd.DataFrame | None = None,
    reference_levels: tuple[str, ...] = (),
) -> str:
    matrix = count_matrix_heatmap.copy()
    matrix.index = matrix.index.map(str)
    matrix.columns = matrix.columns.map(str)
    ordered_sample_ids = [str(sample_id) for sample_id in matrix.columns.tolist()]
    sample_labels, sample_metadata_text = _build_sample_label_map(
        sample_metadata=sample_metadata,
        ordered_sample_ids=ordered_sample_ids,
        reference_levels=reference_levels,
    )

    records = []
    for feature_id in matrix.index.tolist():
        for sample_id in ordered_sample_ids:
            if sample_id not in matrix.columns:
                continue
            records.append(
                {
                    "feature_id": feature_id,
                    "sample_id": sample_id,
                    "sample_label": sample_labels.get(sample_id, sample_id),
                    "sample_metadata": sample_metadata_text.get(sample_id, ""),
                    "value": _value_or_none(matrix.at[feature_id, sample_id]),
                }
            )

    samples = [
        {
            "sample_id": sample_id,
            "sample_label": sample_labels.get(sample_id, sample_id),
            "sample_metadata": sample_metadata_text.get(sample_id, ""),
        }
        for sample_id in ordered_sample_ids
    ]

    payload = {
        "cells": records,
        "feature_order": matrix.index.tolist(),
        "sample_order": ordered_sample_ids,
        "samples": samples,
    }
    return json.dumps(payload).replace("NaN", "null")


def _copy_report_assets(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for folder in ["css", "js", "vega"]:
        src_dir = _ASSETS_DIR / folder
        if src_dir.exists():
            copytree(src_dir, output_dir / folder, dirs_exist_ok=True)


def _render_report(
    output_dir: Path,
    result_table: pd.DataFrame,
    default_effect_id: str,
    default_test_level: str,
    reference_level: str,
    alpha: float,
    include_annotated_results_file: bool,
    sample_distance_matrix: pd.DataFrame | None = None,
    sample_distance_order: tuple[str, ...] = (),
    sample_metadata: pd.DataFrame | None = None,
    reference_levels: tuple[str, ...] = (),
    sample_pca_scores: pd.DataFrame | None = None,
    sample_pca_percent_variance: tuple[float, float] = (),
    count_matrix_heatmap: pd.DataFrame | None = None,
) -> None:
    _copy_report_assets(output_dir)
    data_dir = output_dir / "data"
    data_dir.mkdir(exist_ok=True)
    report_results, default_effect_id = _ensure_effect_columns(
        result_table,
        default_effect_id=default_effect_id,
        default_test_level=default_test_level,
        reference_level=reference_level,
    )
    (data_dir / "results_table.json").write_text(
        _prepare_table_payload(report_results), encoding="utf-8"
    )
    has_sample_distance_heatmap = (
        sample_distance_matrix is not None and not sample_distance_matrix.empty
    )
    has_sample_pca_plot = sample_pca_scores is not None and not sample_pca_scores.empty
    has_count_matrix_heatmap = (
        count_matrix_heatmap is not None and not count_matrix_heatmap.empty
    )
    ordered_sample_ids = []
    if has_sample_distance_heatmap:
        ordered_sample_ids = _resolve_sample_distance_order(
            sample_distance_matrix, sample_distance_order
        )
        (data_dir / "sample_distances.json").write_text(
            _prepare_sample_distance_payload(
                sample_distance_matrix,
                sample_distance_order,
                sample_metadata=sample_metadata,
                reference_levels=reference_levels,
            ),
            encoding="utf-8",
        )
    if has_sample_pca_plot:
        (data_dir / "sample_pca.json").write_text(
            _prepare_sample_pca_payload(
                sample_pca_scores,
                sample_pca_percent_variance=sample_pca_percent_variance,
                sample_distance_order=sample_distance_order,
                sample_metadata=sample_metadata,
                reference_levels=reference_levels,
            ),
            encoding="utf-8",
        )
    if has_count_matrix_heatmap:
        (data_dir / "count_matrix_heatmap.json").write_text(
            _prepare_count_matrix_heatmap_payload(
                count_matrix_heatmap,
                sample_metadata=sample_metadata,
                reference_levels=reference_levels,
            ),
            encoding="utf-8",
        )

    volcano_spec = _load_vega_spec("volcano")
    volcano_spec["signals"][0]["value"] = alpha
    ma_spec = _load_vega_spec("ma")
    ma_spec["signals"][0]["value"] = alpha
    sample_distance_spec = None
    if has_sample_distance_heatmap:
        sample_distance_spec = _load_vega_spec("sample_distance_heatmap")
    sample_pca_spec = None
    if has_sample_pca_plot:
        sample_pca_spec = _load_vega_spec("sample_pca")
    count_matrix_heatmap_spec = None
    if has_count_matrix_heatmap:
        count_matrix_heatmap_spec = _load_vega_spec("count_matrix_heatmap")
    effect_options, default_effect_id, default_effect_label = _collect_effect_options(
        report_results, default_effect_id
    )
    default_results = report_results.loc[
        report_results["effect_id"] == default_effect_id
    ].reset_index(drop=True)
    if default_results.empty:
        default_results = report_results.reset_index(drop=True)

    summary = _summarize_results(default_results, alpha)

    templates = [
        str(_ASSETS_DIR / "index.html"),
    ]
    tabs = [{"title": "Overview", "url": "index.html"}]
    if has_sample_distance_heatmap:
        templates.append(str(_ASSETS_DIR / "sample_distances.html"))
        tabs.append({"title": "Sample distances", "url": "sample_distances.html"})
    templates.append(str(_ASSETS_DIR / "table.html"))
    tabs.append({"title": "Results table", "url": "table.html"})
    context = {
        "tabs": tabs,
        "alpha": str(alpha),
        "summary": summary,
        "vega_volcano_spec": json.dumps(volcano_spec),
        "vega_ma_spec": json.dumps(ma_spec),
        "results_data_path_json": json.dumps("data/results_table.json"),
        "effect_options_json": json.dumps(effect_options),
        "default_effect_id_json": json.dumps(default_effect_id),
        "default_effect_label": default_effect_label,
        "include_annotated_results_file": include_annotated_results_file,
        "has_sample_distance_heatmap": has_sample_distance_heatmap,
        "has_sample_pca_plot": has_sample_pca_plot,
        "has_count_matrix_heatmap": has_count_matrix_heatmap,
        "include_sample_metadata_file": sample_metadata is not None
        and not sample_metadata.empty,
    }
    if has_sample_distance_heatmap:
        context["vega_sample_distance_spec"] = json.dumps(sample_distance_spec)
        context["sample_distances_data_path_json"] = json.dumps(
            "data/sample_distances.json"
        )
        context["sample_distance_order_json"] = json.dumps(ordered_sample_ids)
        context["sample_distance_sample_count"] = len(ordered_sample_ids)
    if has_sample_pca_plot:
        context["vega_sample_pca_spec"] = json.dumps(sample_pca_spec)
        context["sample_pca_data_path_json"] = json.dumps("data/sample_pca.json")
    if has_count_matrix_heatmap:
        context["vega_count_matrix_heatmap_spec"] = json.dumps(
            count_matrix_heatmap_spec
        )
        context["count_matrix_heatmap_data_path_json"] = json.dumps(
            "data/count_matrix_heatmap.json"
        )
        context["count_matrix_heatmap_feature_count"] = len(count_matrix_heatmap.index)
    q2templates.render(templates, str(output_dir), context=context)


def _write_visualization_output(
    output_path: Path,
    run_result: DESeq2RunResult,
    alpha: float,
    display_results: pd.DataFrame | None = None,
    include_annotated_results: bool = False,
) -> None:
    run_result.results.to_csv(output_path / "deseq2_results.tsv", sep="\t", index=False)
    if display_results is None:
        display_results = run_result.results
    elif include_annotated_results:
        display_results.to_csv(
            output_path / "deseq2_results_annotated.tsv", sep="\t", index=False
        )

    run_result.normalized_counts.to_csv(
        output_path / "normalized_counts.tsv", sep="\t", index=False
    )
    if run_result.sample_distance_matrix is not None:
        run_result.sample_distance_matrix.to_csv(
            output_path / "sample_distances.tsv",
            sep="\t",
            index_label="sample_id",
        )
        if (
            run_result.sample_metadata is not None
            and not run_result.sample_metadata.empty
        ):
            sample_metadata = run_result.sample_metadata.copy()
            sample_metadata.index = sample_metadata.index.map(str)
            sample_metadata.columns = sample_metadata.columns.map(str)
            sample_metadata.to_csv(
                output_path / "sample_metadata.tsv",
                sep="\t",
                index_label="sample_id",
            )
    if (
        run_result.count_matrix_heatmap is not None
        and not run_result.count_matrix_heatmap.empty
    ):
        count_matrix_heatmap = run_result.count_matrix_heatmap.copy()
        count_matrix_heatmap.index = count_matrix_heatmap.index.map(str)
        count_matrix_heatmap.columns = count_matrix_heatmap.columns.map(str)
        count_matrix_heatmap.to_csv(
            output_path / "count_matrix_heatmap.tsv",
            sep="\t",
            index_label="feature_id",
        )

    _render_report(
        output_path,
        result_table=display_results,
        default_effect_id=run_result.default_effect_id,
        default_test_level=run_result.test_level,
        reference_level=run_result.reference_level,
        alpha=alpha,
        include_annotated_results_file=include_annotated_results,
        sample_distance_matrix=run_result.sample_distance_matrix,
        sample_distance_order=run_result.sample_distance_order,
        sample_metadata=run_result.sample_metadata,
        reference_levels=run_result.reference_levels,
        sample_pca_scores=run_result.sample_pca_scores,
        sample_pca_percent_variance=run_result.sample_pca_percent_variance,
        count_matrix_heatmap=run_result.count_matrix_heatmap,
    )


def _visualize(
    output_dir: str,
    deseq2_results: DESeq2RunDirectoryFormat,
    gene_annotations: LociDirectoryFormat = None,
    reference_id: str = None,
):
    run_result, alpha = _parse_run_results(deseq2_results)
    annotated_results, include_annotated_results = None, False
    if gene_annotations is not None:
        gff_path = _resolve_reference_gff(gene_annotations, reference_id)
        annotation_table = _load_loci(gff_path)
        annotated_results = _add_annotations(run_result.results, annotation_table)
        include_annotated_results = True

    _write_visualization_output(
        Path(output_dir),
        run_result=run_result,
        alpha=alpha,
        display_results=annotated_results,
        include_annotated_results=include_annotated_results,
    )
