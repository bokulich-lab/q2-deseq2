from unittest.mock import Mock, call

from qiime2.plugin.testing import TestPluginBase

from q2_deseq2 import pipelines


class TestPipelines(TestPluginBase):
    package = "q2_deseq2.tests"

    def test_estimate_without_annotations_calls_actions_with_shared_kwargs(self):
        ctx = Mock()
        table_action = Mock(return_value=("stats-artifact", "run-artifact"))
        visualization_action = Mock(return_value=("viz-artifact",))
        ctx.get_action.side_effect = [table_action, visualization_action]

        observed_stats, observed_visualization = pipelines.estimate(
            ctx=ctx,
            table="table-artifact",
            condition="condition-column",
            gene_annotations=None,
            reference_id=None,
            reference_level="control",
            min_total_count=11,
            fit_type="local",
            alpha=0.01,
            cooks_cutoff=False,
            independent_filtering=False,
        )

        ctx.get_action.assert_has_calls(
            [
                call("deseq2", "_estimate_differential_expression"),
                call("deseq2", "_visualize"),
            ]
        )
        table_action.assert_called_once_with(
            table="table-artifact",
            condition="condition-column",
            reference_level="control",
            min_total_count=11,
            fit_type="local",
            alpha=0.01,
            cooks_cutoff=False,
            independent_filtering=False,
        )
        visualization_action.assert_called_once_with(deseq2_results="run-artifact")
        self.assertEqual(observed_stats, "stats-artifact")
        self.assertEqual(observed_visualization, "viz-artifact")

    def test_estimate_with_annotations_passes_annotation_inputs_to_visualizer(self):
        ctx = Mock()
        table_action = Mock(return_value=("stats-artifact", "run-artifact"))
        visualization_action = Mock(return_value=("viz-artifact",))
        ctx.get_action.side_effect = [table_action, visualization_action]

        observed_stats, observed_visualization = pipelines.estimate(
            ctx=ctx,
            table="table-artifact",
            condition="condition-column",
            gene_annotations="gff-artifact",
            reference_id="ref-a",
        )

        table_action.assert_called_once_with(
            table="table-artifact",
            condition="condition-column",
            reference_level="",
            min_total_count=10,
            fit_type="parametric",
            alpha=0.05,
            cooks_cutoff=True,
            independent_filtering=True,
        )
        visualization_action.assert_called_once_with(
            deseq2_results="run-artifact",
            gene_annotations="gff-artifact",
            reference_id="ref-a",
        )
        self.assertEqual(observed_stats, "stats-artifact")
        self.assertEqual(observed_visualization, "viz-artifact")
