import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_deseq2._run_data import DESeq2RunResult
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
        self.run_result = DESeq2RunResult(
            results=self.multi_results,
            normalized_counts=self.normalized_counts,
            ma_plot_png=b"ma-plot",
            volcano_plot_png=b"volcano-plot",
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
            self.assertTrue((output_path / "data" / "results_table.json").exists())
            self.assertTrue((output_path / "css" / "styles.css").exists())
            self.assertTrue((output_path / "js" / "linked_plots.js").exists())
            self.assertTrue((output_path / "vega" / "volcano.json").exists())
            self.assertTrue((output_path / "vega" / "ma.json").exists())

            report_payload = json.loads(
                (output_path / "data" / "results_table.json").read_text(
                    encoding="utf-8"
                )
            )

        self.assertEqual(captured["output_dir"], temp_dir)
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
        self.assertEqual(
            captured["context"]["default_effect_label"], "condition: other vs control"
        )
        self.assertTrue(captured["context"]["include_annotated_results_file"])
        self.assertEqual(captured["context"]["summary"]["total_features"], 3)
        self.assertEqual(captured["context"]["summary"]["significant_features"], 1)
        self.assertEqual(captured["context"]["summary"]["annotated_features"], 3)
        self.assertEqual(
            len(json.loads(captured["context"]["effect_options_json"])), 2
        )
        self.assertEqual(report_payload["columns"][0], "feature_id")
        self.assertEqual(len(report_payload["data"]), 6)
        self.assertIn("effect_id", report_payload["columns"])

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
