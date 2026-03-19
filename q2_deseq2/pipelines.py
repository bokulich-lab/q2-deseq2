# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
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
    test_level="",
    reference_level="",
    min_total_count=10,
    fit_type="parametric",
    alpha=0.05,
    cooks_cutoff=True,
    independent_filtering=True,
):
    table_action = ctx.get_action("deseq2", "_differential_expression_table")
    visualization_action = ctx.get_action("deseq2", "_differential_expression")

    shared_kwargs = {
        "table": table,
        "condition": condition,
        "test_level": test_level,
        "reference_level": reference_level,
        "min_total_count": min_total_count,
        "fit_type": fit_type,
        "alpha": alpha,
        "cooks_cutoff": cooks_cutoff,
        "independent_filtering": independent_filtering,
    }

    differential_expression_stats, differential_expression_run = table_action(
        **shared_kwargs
    )
    visualization_kwargs = {"differential_expression_run": differential_expression_run}
    if annotations is None:
        (differential_expression_visualization,) = visualization_action(
            **visualization_kwargs
        )
    else:
        visualization_kwargs["annotations"] = annotations
        (differential_expression_visualization,) = visualization_action(
            **visualization_kwargs
        )

    return differential_expression_stats, differential_expression_visualization
