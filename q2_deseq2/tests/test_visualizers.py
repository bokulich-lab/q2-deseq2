import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd
from q2_deseq2._run_data import DESeq2RunResult
from qiime2.plugin.testing import TestPluginBase

from q2_deseq2 import visualizers


class _FakeAnnotations:
    def __init__(self, file_map):
        self._file_map = file_map

    def file_dict(self):
        return self._file_map


class TestVisualizers(TestPluginBase):
    package = "q2_deseq2.tests"

    def setUp(self):
        super().setUp()
        self.raw_results = pd.DataFrame(
            [
                {
                    "feature_id": "GG_OTU_1",
                    "baseMean": 120.0,
                    "log2FoldChange": 2.4,
                    "lfcSE": 0.3,
                    "stat": 4.7,
                    "pvalue": 0.002,
                    "padj": 0.01,
                },
                {
                    "feature_id": "GG_OTU_2",
                    "baseMean": 55.0,
                    "log2FoldChange": -0.8,
                    "lfcSE": 0.2,
                    "stat": -3.4,
                    "pvalue": 0.3,
                    "padj": 0.5,
                },
                {
                    "feature_id": "GG_OTU_3",
                    "baseMean": 10.0,
                    "log2FoldChange": None,
                    "lfcSE": 0.1,
                    "stat": None,
                    "pvalue": None,
                    "padj": None,
                },
            ]
        )
        self.annotated_results = self.raw_results.assign(
            gene_name=["dnaA", "gene-b", None],
            product=["Replication initiator", "Some product", "Third product"],
        )
        self.multi_results = pd.concat(
            [
                self.raw_results.assign(
                    effect_id="contrast::condition::other::control",
                    effect_label="condition: other vs control",
                    effect_kind="contrast",
                    effect_expression='contrast=c("condition","other","control")',
                    comparison="other vs. control",
                    test_level="other",
                    reference_level="control",
                ),
                self.raw_results.assign(
                    effect_id="contrast::condition::treated::control",
                    effect_label="condition: treated vs control",
                    effect_kind="contrast",
                    effect_expression='contrast=c("condition","treated","control")',
                    comparison="treated vs. control",
                    test_level="treated",
                    reference_level="control",
                ),
            ],
            ignore_index=True,
        )
        self.annotated_multi_results = pd.concat(
            [
                self.annotated_results.assign(
                    effect_id="contrast::condition::other::control",
                    effect_label="condition: other vs control",
                    effect_kind="contrast",
                    effect_expression='contrast=c("condition","other","control")',
                    comparison="other vs. control",
                    test_level="other",
                    reference_level="control",
                ),
                self.annotated_results.assign(
                    effect_id="contrast::condition::treated::control",
                    effect_label="condition: treated vs control",
                    effect_kind="contrast",
                    effect_expression='contrast=c("condition","treated","control")',
                    comparison="treated vs. control",
                    test_level="treated",
                    reference_level="control",
                ),
            ],
            ignore_index=True,
        )
        self.normalized_counts = pd.DataFrame(
            [
                {"feature_id": "GG_OTU_1", "Sample1": 100.0, "Sample2": 140.0},
                {"feature_id": "GG_OTU_2", "Sample1": 50.0, "Sample2": 60.0},
            ]
        )
        self.sample_distance_matrix = pd.DataFrame(
            [[0.0, 1.75], [1.75, 0.0]],
            index=["Sample1", "Sample2"],
            columns=["Sample1", "Sample2"],
        )
        self.sample_metadata = pd.DataFrame(
            {
                "batch": ["A", "B"],
                "condition": ["control", "treated"],
                "depth": [11.0, 14.5],
            },
            index=["Sample1", "Sample2"],
        )
        self.run_result = DESeq2RunResult(
            results=self.multi_results,
            normalized_counts=self.normalized_counts,
            sample_metadata=self.sample_metadata,
            test_level="other",
            reference_level="control",
            default_effect_id="contrast::condition::other::control",
            fixed_effects_formula="condition",
            reference_levels=("condition::control",),
            test="wald",
            reduced_formula="",
            available_results_names=("Intercept", "condition_treated_vs_control"),
            selected_effect_specs=(
                "contrast::condition::other::control",
                "contrast::condition::treated::control",
            ),
            sample_distance_matrix=self.sample_distance_matrix,
            sample_distance_order=("Sample2", "Sample1"),
            sample_pca_scores=pd.DataFrame(
                {"PC1": [-2.1, 2.1], "PC2": [0.35, -0.35]},
                index=["Sample1", "Sample2"],
            ),
            sample_pca_percent_variance=(68.2, 21.5),
            count_matrix_heatmap=pd.DataFrame(
                {"Sample2": [7.5, 5.1], "Sample1": [6.8, 4.6]},
                index=["GG_OTU_2", "GG_OTU_1"],
            ),
        )

    def test_value_or_none_handles_missing_and_numeric_values(self):
        self.assertIsNone(visualizers._value_or_none(pd.NA))
        self.assertEqual(visualizers._value_or_none("2.5"), 2.5)
        self.assertIsNone(visualizers._value_or_none("not-a-number"))

    def test_parse_gff3_attributes_and_first_non_empty(self):
        observed = visualizers._parse_gff3_attributes(
            "ID=feat1;Name=dnaA;product=Replication%20initiator;flag"
        )

        self.assertEqual(
            observed,
            {
                "ID": "feat1",
                "Name": "dnaA",
                "product": "Replication initiator",
                "flag": "",
            },
        )
        self.assertEqual(
            visualizers._first_non_empty(
                {"gene": "", "Name": "abc,def", "product": "xyz"},
                ["gene", "Name", "product"],
            ),
            "abc",
        )
        self.assertIsNone(visualizers._first_non_empty({}, ["gene"]))

    def test_resolve_reference_gff_returns_matching_path(self):
        gff_path = Path(self.get_data_path("reference.gff3"))
        annotations = _FakeAnnotations({"ref-a": gff_path})

        observed = visualizers._resolve_reference_gff(annotations, "ref-a")

        self.assertEqual(observed, gff_path)

    def test_resolve_reference_gff_raises_for_missing_reference(self):
        annotations = _FakeAnnotations(
            {"ref-a": Path(self.get_data_path("reference.gff3"))}
        )

        with self.assertRaisesRegex(ValueError, 'Reference "missing" was not found'):
            visualizers._resolve_reference_gff(annotations, "missing")

    def test_load_loci_parses_reference_fixture(self):
        observed = visualizers._load_loci(Path(self.get_data_path("reference.gff3")))
        observed = observed.set_index("feature_id")

        self.assertEqual(observed.loc["GG_OTU_1", "gene_name"], "dnaA")
        self.assertEqual(observed.loc["GG_OTU_1", "product"], "Replication initiator")
        self.assertEqual(observed.loc["GG_OTU_2", "gene_name"], "gene-b")
        self.assertEqual(observed.loc["GG_OTU_2", "product"], "Some product")
        self.assertEqual(observed.loc["GG_OTU_3", "gene_name"], "GG_OTU_3")
        self.assertEqual(observed.loc["GG_OTU_3", "product"], "Third product")
        self.assertIn("dnaA", observed.index)
        self.assertIn("gene-b", observed.index)

    def test_add_annotations_merges_on_feature_id(self):
        annotation_table = pd.DataFrame(
            [
                {
                    "feature_id": "GG_OTU_1",
                    "gene_name": "dnaA",
                    "product": "Replication",
                },
                {"feature_id": "GG_OTU_9", "gene_name": "unused", "product": "unused"},
            ]
        )

        observed = visualizers._add_annotations(self.raw_results, annotation_table)

        self.assertEqual(observed.loc[0, "gene_name"], "dnaA")
        self.assertTrue(pd.isna(observed.loc[1, "gene_name"]))
        self.assertTrue(pd.isna(observed.loc[2, "gene_name"]))
        self.assertEqual(observed.loc[0, "product"], "Replication")
        self.assertTrue(pd.isna(observed.loc[1, "product"]))

    def test_build_plot_records_serializes_nullable_fields(self):
        observed = visualizers._build_plot_records(self.annotated_results)

        self.assertEqual(observed[0]["feature_id"], "GG_OTU_1")
        self.assertEqual(observed[0]["gene_name"], "dnaA")
        self.assertEqual(observed[0]["product"], "Replication initiator")
        self.assertEqual(observed[0]["log2FoldChange"], 2.4)
        self.assertIsNone(observed[2]["log2FoldChange"])
        self.assertIsNone(observed[2]["gene_name"])

    def test_ensure_effect_columns_backfills_legacy_single_comparison_results(self):
        observed, default_effect_id = visualizers._ensure_effect_columns(
            self.annotated_results,
            default_effect_id="",
            default_test_level="treated",
            reference_level="control",
        )

        self.assertTrue(
            {
                "effect_id",
                "effect_label",
                "effect_kind",
                "effect_expression",
                "comparison",
                "test_level",
                "reference_level",
            }
            <= set(observed.columns)
        )
        self.assertEqual(observed["effect_id"].nunique(), 1)
        self.assertEqual(observed.loc[0, "effect_label"], "treated vs. control")
        self.assertEqual(observed.loc[0, "effect_kind"], "comparison")
        self.assertEqual(default_effect_id, "legacy::treated::control")

    def test_collect_effect_options_returns_default_label(self):
        effect_options, default_effect_id, default_effect_label = (
            visualizers._collect_effect_options(
                self.multi_results, "contrast::condition::other::control"
            )
        )

        self.assertEqual(default_effect_id, "contrast::condition::other::control")
        self.assertEqual(default_effect_label, "condition: other vs control")
        self.assertEqual(len(effect_options), 2)

    def test_summarize_results_counts_expected_categories(self):
        observed = visualizers._summarize_results(self.annotated_results, alpha=0.05)

        self.assertEqual(
            observed,
            {
                "total_features": 3,
                "plottable_features": 2,
                "significant_features": 1,
                "alpha_display": "0.05",
                "annotated_features": 3,
            },
        )

    def test_prepare_table_payload_serializes_nulls(self):
        payload = json.loads(visualizers._prepare_table_payload(self.annotated_results))

        self.assertEqual(payload["columns"][0], "feature_id")
        self.assertEqual(payload["data"][0][0], "GG_OTU_1")
        self.assertIsNone(payload["data"][2][2])

    def test_prepare_sample_distance_payload_preserves_cluster_order(self):
        payload = json.loads(
            visualizers._prepare_sample_distance_payload(
                self.sample_distance_matrix,
                ("Sample2", "Sample1"),
                sample_metadata=self.sample_metadata,
                reference_levels=("condition::control",),
            )
        )

        self.assertEqual(payload[0]["sample_x"], "Sample2")
        self.assertEqual(payload[0]["sample_y"], "Sample2")
        self.assertEqual(
            payload[0]["sample_y_label"],
            "Sample2 | condition=treated | batch=B",
        )
        self.assertEqual(payload[1]["sample_x"], "Sample1")
        self.assertEqual(payload[1]["sample_y"], "Sample2")
        self.assertEqual(payload[1]["distance"], 1.75)
        self.assertEqual(payload[1]["sample_x_metadata"], "condition=control; batch=A")
        self.assertEqual(payload[1]["sample_y_metadata"], "condition=treated; batch=B")

    def test_build_sample_label_map_prioritizes_reference_columns(self):
        sample_labels, sample_metadata_text = visualizers._build_sample_label_map(
            self.sample_metadata,
            ["Sample2", "Sample1"],
            ("condition::control",),
        )

        self.assertEqual(
            sample_labels["Sample2"],
            "Sample2 | condition=treated | batch=B",
        )
        self.assertEqual(
            sample_metadata_text["Sample1"],
            "condition=control; batch=A",
        )

    def test_prepare_sample_pca_payload_uses_reference_grouping(self):
        payload = json.loads(
            visualizers._prepare_sample_pca_payload(
                self.run_result.sample_pca_scores,
                sample_pca_percent_variance=self.run_result.sample_pca_percent_variance,
                sample_distance_order=("Sample2", "Sample1"),
                sample_metadata=self.sample_metadata,
                reference_levels=("condition::control",),
            )
        )

        self.assertEqual(payload["group_field"], "condition")
        self.assertEqual(payload["group_label"], "condition (ref: control)")
        self.assertEqual(payload["percent_variance"]["PC1"], 68.2)
        self.assertEqual(payload["points"][0]["sample_id"], "Sample2")
        self.assertEqual(payload["points"][0]["group_value"], "treated")
        self.assertEqual(
            payload["points"][1]["sample_metadata"],
            "condition=control; batch=A",
        )

    def test_prepare_count_matrix_heatmap_payload_preserves_order(self):
        payload = json.loads(
            visualizers._prepare_count_matrix_heatmap_payload(
                self.run_result.count_matrix_heatmap,
                sample_metadata=self.sample_metadata,
                reference_levels=("condition::control",),
            )
        )

        self.assertEqual(payload["feature_order"], ["GG_OTU_2", "GG_OTU_1"])
        self.assertEqual(payload["sample_order"], ["Sample2", "Sample1"])
        self.assertEqual(payload["samples"][0]["sample_id"], "Sample2")
        self.assertEqual(
            payload["samples"][0]["sample_metadata"],
            "condition=treated; batch=B",
        )
        self.assertEqual(payload["cells"][0]["feature_id"], "GG_OTU_2")
        self.assertEqual(payload["cells"][0]["sample_id"], "Sample2")
        self.assertEqual(payload["cells"][0]["value"], 7.5)
        self.assertEqual(
            payload["cells"][0]["sample_metadata"], "condition=treated; batch=B"
        )

    def test_write_visualization_output_renders_tabbed_report(self):
        captured = {}

        def fake_render(templates, output_dir, context):
            captured["templates"] = templates
            captured["output_dir"] = output_dir
            captured["context"] = context

        with TemporaryDirectory() as temp_dir, patch.object(
            visualizers, "q2templates", SimpleNamespace(render=fake_render)
        ):
            output_path = Path(temp_dir)
            visualizers._write_visualization_output(
                output_path,
                run_result=self.run_result,
                alpha=0.05,
                display_results=self.annotated_multi_results,
                include_annotated_results=True,
            )

            self.assertTrue((output_path / "deseq2_results.tsv").exists())
            self.assertTrue((output_path / "deseq2_results_annotated.tsv").exists())
            self.assertTrue((output_path / "normalized_counts.tsv").exists())
            self.assertTrue((output_path / "sample_distances.tsv").exists())
            self.assertTrue((output_path / "sample_metadata.tsv").exists())
            self.assertTrue((output_path / "count_matrix_heatmap.tsv").exists())
            self.assertFalse((output_path / "ma_plot.png").exists())
            self.assertFalse((output_path / "volcano_plot.png").exists())
            self.assertTrue((output_path / "data" / "results_table.json").exists())
            self.assertTrue((output_path / "data" / "sample_distances.json").exists())
            self.assertTrue((output_path / "data" / "sample_pca.json").exists())
            self.assertTrue(
                (output_path / "data" / "count_matrix_heatmap.json").exists()
            )
            self.assertTrue((output_path / "css" / "styles.css").exists())
            self.assertTrue((output_path / "js" / "linked_plots.js").exists())
            self.assertTrue((output_path / "js" / "sample_distances.js").exists())
            self.assertTrue((output_path / "vega" / "volcano.json").exists())
            self.assertTrue((output_path / "vega" / "ma.json").exists())
            self.assertTrue(
                (output_path / "vega" / "sample_distance_heatmap.json").exists()
            )
            self.assertTrue((output_path / "vega" / "sample_pca.json").exists())
            self.assertTrue(
                (output_path / "vega" / "count_matrix_heatmap.json").exists()
            )

            report_payload = json.loads(
                (output_path / "data" / "results_table.json").read_text(
                    encoding="utf-8"
                )
            )
            sample_distance_payload = json.loads(
                (output_path / "data" / "sample_distances.json").read_text(
                    encoding="utf-8"
                )
            )
            sample_pca_payload = json.loads(
                (output_path / "data" / "sample_pca.json").read_text(encoding="utf-8")
            )
            count_matrix_heatmap_payload = json.loads(
                (output_path / "data" / "count_matrix_heatmap.json").read_text(
                    encoding="utf-8"
                )
            )

        self.assertEqual(captured["output_dir"], temp_dir)
        self.assertEqual(
            [Path(template).name for template in captured["templates"]],
            ["index.html", "sample_distances.html", "table.html"],
        )
        self.assertEqual(
            captured["context"]["tabs"],
            [
                {"title": "Overview", "url": "index.html"},
                {"title": "Sample distances", "url": "sample_distances.html"},
                {"title": "Results table", "url": "table.html"},
            ],
        )
        self.assertEqual(
            captured["context"]["default_effect_label"], "condition: other vs control"
        )
        self.assertTrue(captured["context"]["include_annotated_results_file"])
        self.assertEqual(captured["context"]["summary"]["total_features"], 3)
        self.assertEqual(captured["context"]["summary"]["significant_features"], 1)
        self.assertEqual(captured["context"]["summary"]["annotated_features"], 3)
        self.assertEqual(len(json.loads(captured["context"]["effect_options_json"])), 2)
        self.assertTrue(captured["context"]["include_sample_metadata_file"])
        self.assertTrue(captured["context"]["has_sample_pca_plot"])
        self.assertTrue(captured["context"]["has_count_matrix_heatmap"])
        self.assertEqual(report_payload["columns"][0], "feature_id")
        self.assertEqual(len(report_payload["data"]), 6)
        self.assertIn("effect_id", report_payload["columns"])
        self.assertEqual(len(sample_distance_payload), 4)
        self.assertEqual(
            sample_distance_payload[0]["sample_y_label"],
            "Sample2 | condition=treated | batch=B",
        )
        self.assertEqual(
            json.loads(captured["context"]["sample_distance_order_json"]),
            ["Sample2", "Sample1"],
        )
        self.assertEqual(captured["context"]["sample_distance_sample_count"], 2)
        self.assertEqual(sample_pca_payload["group_label"], "condition (ref: control)")
        self.assertEqual(sample_pca_payload["points"][0]["sample_id"], "Sample2")
        self.assertEqual(
            sample_pca_payload["percent_variance"],
            {"PC1": 68.2, "PC2": 21.5},
        )
        self.assertEqual(
            json.loads(captured["context"]["sample_pca_data_path_json"]),
            "data/sample_pca.json",
        )
        self.assertEqual(
            json.loads(captured["context"]["count_matrix_heatmap_data_path_json"]),
            "data/count_matrix_heatmap.json",
        )
        self.assertEqual(
            count_matrix_heatmap_payload["feature_order"],
            ["GG_OTU_2", "GG_OTU_1"],
        )

    def test_write_visualization_output_skips_sample_distance_tab_without_matrix(self):
        captured = {}

        def fake_render(templates, output_dir, context):
            captured["templates"] = templates
            captured["output_dir"] = output_dir
            captured["context"] = context

        run_result = self.run_result._replace(
            sample_distance_matrix=None,
            sample_distance_order=(),
            sample_pca_scores=None,
            sample_pca_percent_variance=(),
            count_matrix_heatmap=None,
        )

        with TemporaryDirectory() as temp_dir, patch.object(
            visualizers, "q2templates", SimpleNamespace(render=fake_render)
        ):
            output_path = Path(temp_dir)
            visualizers._write_visualization_output(
                output_path,
                run_result=run_result,
                alpha=0.05,
                display_results=self.annotated_multi_results,
                include_annotated_results=True,
            )

            self.assertFalse((output_path / "sample_distances.tsv").exists())
            self.assertFalse((output_path / "sample_metadata.tsv").exists())
            self.assertFalse((output_path / "count_matrix_heatmap.tsv").exists())
            self.assertFalse((output_path / "data" / "sample_distances.json").exists())
            self.assertFalse((output_path / "data" / "sample_pca.json").exists())
            self.assertFalse(
                (output_path / "data" / "count_matrix_heatmap.json").exists()
            )

        self.assertEqual(
            [Path(template).name for template in captured["templates"]],
            ["index.html", "table.html"],
        )
        self.assertEqual(
            captured["context"]["tabs"],
            [
                {"title": "Overview", "url": "index.html"},
                {"title": "Results table", "url": "table.html"},
            ],
        )

    @patch("q2_deseq2.visualizers._write_visualization_output")
    @patch("q2_deseq2.visualizers._parse_run_results")
    def test_visualize_without_annotations_calls_writer_with_raw_results(
        self, parse_run_results_mock, write_output_mock
    ):
        parse_run_results_mock.return_value = (self.run_result, 0.05)

        visualizers._visualize(
            output_dir="fake-output",
            deseq2_results="run-data",
            gene_annotations=None,
            reference_id=None,
        )

        parse_run_results_mock.assert_called_once_with("run-data")
        write_output_mock.assert_called_once_with(
            Path("fake-output"),
            run_result=self.run_result,
            alpha=0.05,
            display_results=None,
            include_annotated_results=False,
        )

    @patch("q2_deseq2.visualizers._write_visualization_output")
    @patch("q2_deseq2.visualizers._add_annotations")
    @patch("q2_deseq2.visualizers._load_loci")
    @patch("q2_deseq2.visualizers._resolve_reference_gff")
    @patch("q2_deseq2.visualizers._parse_run_results")
    def test_visualize_with_annotations_merges_and_writes(
        self,
        parse_run_results_mock,
        resolve_reference_gff_mock,
        load_loci_mock,
        add_annotations_mock,
        write_output_mock,
    ):
        parse_run_results_mock.return_value = (self.run_result, 0.05)
        resolve_reference_gff_mock.return_value = Path("/tmp/reference.gff3")
        load_loci_mock.return_value = pd.DataFrame(
            [{"feature_id": "GG_OTU_1", "gene_name": "dnaA", "product": "Replication"}]
        )
        add_annotations_mock.return_value = self.annotated_multi_results

        visualizers._visualize(
            output_dir="fake-output",
            deseq2_results="run-data",
            gene_annotations="annotations",
            reference_id="ref-a",
        )

        parse_run_results_mock.assert_called_once_with("run-data")
        resolve_reference_gff_mock.assert_called_once_with("annotations", "ref-a")
        load_loci_mock.assert_called_once_with(Path("/tmp/reference.gff3"))
        add_annotations_mock.assert_called_once()
        write_output_mock.assert_called_once_with(
            Path("fake-output"),
            run_result=self.run_result,
            alpha=0.05,
            display_results=self.annotated_multi_results,
            include_annotated_results=True,
        )
