# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def estimate(
    ctx,
    table,
    condition,
    gene_annotations=None,
    reference_id=None,
    reference_level="",
    min_total_count=10,
    fit_type="parametric",
    size_factor_type="poscounts",
    alpha=0.05,
    cooks_cutoff=True,
    independent_filtering=True,
):
    table_action = ctx.get_action("deseq2", "_estimate_differential_expression")
    visualization_action = ctx.get_action("deseq2", "_visualize")

    shared_kwargs = {
        "table": table,
        "condition": condition,
        "reference_level": reference_level,
        "min_total_count": min_total_count,
        "fit_type": fit_type,
        "size_factor_type": size_factor_type,
        "alpha": alpha,
        "cooks_cutoff": cooks_cutoff,
        "independent_filtering": independent_filtering,
    }

    differential_expression_stats, differential_expression_run = table_action(
        **shared_kwargs
    )
    visualization_kwargs = {"deseq2_results": differential_expression_run}
    if gene_annotations is None:
        (differential_expression_visualization,) = visualization_action(
            **visualization_kwargs
        )
    else:
        visualization_kwargs["gene_annotations"] = gene_annotations
        visualization_kwargs["reference_id"] = reference_id
        (differential_expression_visualization,) = visualization_action(
            **visualization_kwargs
        )

    return differential_expression_stats, differential_expression_visualization


def estimate_model(
    ctx,
    table,
    metadata,
    fixed_effects_formula,
    gene_annotations=None,
    reference_id=None,
    reference_levels=None,
    effect_specs=None,
    test="wald",
    reduced_formula="",
    min_total_count=10,
    fit_type="parametric",
    size_factor_type="poscounts",
    alpha=0.05,
    cooks_cutoff=True,
    independent_filtering=True,
):
    table_action = ctx.get_action("deseq2", "_estimate_model")
    visualization_action = ctx.get_action("deseq2", "_visualize")

    shared_kwargs = {
        "table": table,
        "metadata": metadata,
        "fixed_effects_formula": fixed_effects_formula,
        "reference_levels": reference_levels or [""],
        "effect_specs": effect_specs or [""],
        "test": test,
        "reduced_formula": reduced_formula,
        "min_total_count": min_total_count,
        "fit_type": fit_type,
        "size_factor_type": size_factor_type,
        "alpha": alpha,
        "cooks_cutoff": cooks_cutoff,
        "independent_filtering": independent_filtering,
    }

    differential_expression_stats, differential_expression_run = table_action(
        **shared_kwargs
    )
    visualization_kwargs = {"deseq2_results": differential_expression_run}
    if gene_annotations is None:
        (differential_expression_visualization,) = visualization_action(
            **visualization_kwargs
        )
    else:
        visualization_kwargs["gene_annotations"] = gene_annotations
        visualization_kwargs["reference_id"] = reference_id
        (differential_expression_visualization,) = visualization_action(
            **visualization_kwargs
        )

    return differential_expression_stats, differential_expression_visualization
