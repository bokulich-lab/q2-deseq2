# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import pandas as pd
from q2_deseq2._frame_utils import (
    _first_non_empty_string,
    _first_value_from_column,
    _unique_non_empty_values,
)
from q2_deseq2._run_data import DESeq2RunResult, write_run_result_artifact

from q2_deseq2.types import DESeq2RunDirectoryFormat
from q2_deseq2.utils.analytics import (
    _compute_count_matrix_heatmap,
    _compute_normalized_counts,
    _compute_run_analytics,
    _compute_sample_distance_matrix,
    _compute_sample_distance_order,
    _compute_sample_pca,
)
from q2_deseq2.utils.prep import (
    _coerce_metadata_column,
    _coerce_metadata_frame,
    _collect_matching_samples,
    _extract_formula_columns,
    _filter_counts,
    _normalize_formula,
    _normalize_reference_levels,
    _normalize_size_factor_type,
    _parse_reference_level_spec,
    _prepare_count_table,
    _prepare_inputs,
    _prepare_model_inputs,
    _validate_effect_specs,
)
from q2_deseq2.utils.runner import _run_deseq2_with_frames, _write_lines

__all__ = [
    "_coerce_metadata_column",
    "_coerce_metadata_frame",
    "_collect_matching_samples",
    "_compute_count_matrix_heatmap",
    "_compute_normalized_counts",
    "_compute_run_analytics",
    "_compute_sample_distance_matrix",
    "_compute_sample_distance_order",
    "_compute_sample_pca",
    "_estimate_differential_expression",
    "_estimate_model",
    "_extract_formula_columns",
    "_filter_counts",
    "_first_non_empty_string",
    "_first_value_from_column",
    "_normalize_formula",
    "_normalize_reference_levels",
    "_normalize_size_factor_type",
    "_parse_reference_level_spec",
    "_prepare_count_table",
    "_prepare_inputs",
    "_prepare_model_inputs",
    "_run_deseq2_with_frames",
    "_unique_non_empty_values",
    "_validate_effect_specs",
    "_write_lines",
    "DESeq2RunResult",
    "run_deseq2",
    "run_deseq2_model",
    "write_run_result_artifact",
]


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
    _, coldata_df, comparison_levels, reference_level = _prepare_inputs(
        table, condition, min_total_count, reference_level
    )
    effect_specs = [
        f"contrast::condition::{comparison_level}::{reference_level}"
        for comparison_level in comparison_levels
    ]

    run_result = run_deseq2_model(
        table=table,
        metadata=coldata_df,
        fixed_effects_formula="condition",
        reference_levels=[f"condition::{reference_level}"],
        effect_specs=effect_specs,
        test="wald",
        reduced_formula="",
        min_total_count=min_total_count,
        fit_type=fit_type,
        size_factor_type=size_factor_type,
        alpha=alpha,
        cooks_cutoff=cooks_cutoff,
        independent_filtering=independent_filtering,
    )
    return run_result._replace(
        test_level=comparison_levels[0] if len(comparison_levels) == 1 else "",
        reference_level=reference_level,
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
