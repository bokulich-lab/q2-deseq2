# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile
from pathlib import Path
from subprocess import CalledProcessError, run

import pandas as pd
from q2_deseq2._frame_utils import _first_value_from_column, _unique_non_empty_values
from q2_deseq2._run_data import DESeq2RunResult

from q2_deseq2.utils.analytics import _compute_run_analytics
from q2_deseq2.utils.prep import (
    _normalize_formula,
    _normalize_reference_levels,
    _normalize_size_factor_type,
    _validate_effect_specs,
)


def _write_lines(path: Path, lines: list[str]) -> None:
    payload = "\n".join(lines)
    if payload:
        payload += "\n"
    path.write_text(payload, encoding="utf-8")


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
        script_fp = Path(__file__).resolve().parent.parent / "r" / "run_deseq2.R"
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
