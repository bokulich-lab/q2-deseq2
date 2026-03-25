from pathlib import Path
from subprocess import CalledProcessError
from tempfile import TemporaryDirectory
from unittest.mock import Mock, patch

import biom
import pandas as pd
import qiime2
from pandas.testing import assert_frame_equal
from qiime2.plugin.testing import TestPluginBase

from q2_deseq2._run_data import DESeq2RunResult
from q2_deseq2 import methods


class TestMethods(TestPluginBase):
    package = "q2_deseq2.tests"

    def setUp(self):
        super().setUp()
        self.table = biom.load_table(self.get_data_path("table-1.biom"))
        self.condition = qiime2.Metadata.load(
            self.get_data_path("condition.tsv")
        ).get_column("condition")
        self.condition_single_level = qiime2.Metadata.load(
            self.get_data_path("condition-single-level.tsv")
        ).get_column("condition")
        self.condition_three_levels = qiime2.Metadata.load(
            self.get_data_path("condition-three-levels.tsv")
        ).get_column("condition")
        self.condition_unmatched = qiime2.Metadata.load(
            self.get_data_path("condition-unmatched.tsv")
        ).get_column("condition")
        self.model_metadata = qiime2.Metadata.load(
            self.get_data_path("model-metadata.tsv")
        )
        self.sample_run_result = DESeq2RunResult(
            results=pd.DataFrame(
                [
                    {
                        "effect_id": "contrast::condition::other::control",
                        "effect_label": "condition: other vs control",
                        "effect_kind": "contrast",
                        "effect_expression": 'contrast=c("condition","other","control")',
                        "comparison": "other vs. control",
                        "test_level": "other",
                        "reference_level": "control",
                        "feature_id": "GG_OTU_2",
                        "baseMean": 120.0,
                        "log2FoldChange": 2.4,
                        "lfcSE": 0.4,
                        "stat": 6.0,
                        "pvalue": 0.002,
                        "padj": 0.01,
                    },
                    {
                        "effect_id": "contrast::condition::treated::control",
                        "effect_label": "condition: treated vs control",
                        "effect_kind": "contrast",
                        "effect_expression": 'contrast=c("condition","treated","control")',
                        "comparison": "treated vs. control",
                        "test_level": "treated",
                        "reference_level": "control",
                        "feature_id": "GG_OTU_2",
                        "baseMean": 140.0,
                        "log2FoldChange": 1.2,
                        "lfcSE": 0.3,
                        "stat": 4.0,
                        "pvalue": 0.01,
                        "padj": 0.04,
                    },
                ]
            ),
            normalized_counts=pd.DataFrame(
                [{"feature_id": "GG_OTU_2", "Sample1": 3.0, "Sample2": 4.0}]
            ),
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
            sample_distance_matrix=pd.DataFrame(
                [[0.0, 1.25], [1.25, 0.0]],
                index=["Sample1", "Sample2"],
                columns=["Sample1", "Sample2"],
            ),
            sample_distance_order=("Sample2", "Sample1"),
        )

    def test_prepare_inputs_infers_reference_for_two_level_metadata(self):
        observed_counts, observed_coldata, observed_comparisons, observed_reference = (
            methods._prepare_inputs(
                self.table,
                self.condition,
                min_total_count=3,
                reference_level="",
            )
        )

        expected_counts = self.table.to_dataframe(dense=True).round().astype(int)
        expected_counts = expected_counts.loc[expected_counts.sum(axis=1) >= 3, :]
        expected_coldata = pd.DataFrame(
            {"condition": self.condition.to_series().astype(str)}
        )

        assert_frame_equal(observed_counts, expected_counts)
        assert_frame_equal(observed_coldata, expected_coldata)
        self.assertEqual(observed_comparisons, ["treated"])
        self.assertEqual(observed_reference, "control")
        self.assertTrue(
            all(dtype.kind in {"i", "u"} for dtype in observed_counts.dtypes)
        )

    def test_prepare_inputs_accepts_reference_level_for_three_level_metadata(self):
        observed_counts, observed_coldata, observed_comparisons, observed_reference = (
            methods._prepare_inputs(
                self.table,
                self.condition_three_levels,
                min_total_count=0,
                reference_level="other",
            )
        )

        expected_counts = self.table.to_dataframe(dense=True).round().astype(int)
        expected_coldata = pd.DataFrame(
            {"condition": self.condition_three_levels.to_series().astype(str)}
        )

        assert_frame_equal(observed_counts, expected_counts)
        assert_frame_equal(observed_coldata, expected_coldata)
        self.assertEqual(observed_comparisons, ["control", "treated"])
        self.assertEqual(observed_reference, "other")

    def test_prepare_model_inputs_accepts_formula_and_reference_levels(self):
        (
            observed_counts,
            observed_coldata,
            observed_formula,
            observed_reference_levels,
            observed_reduced,
        ) = methods._prepare_model_inputs(
            self.table,
            self.model_metadata,
            fixed_effects_formula="genotype + treatment + genotype:treatment",
            min_total_count=0,
            reference_levels=["genotype::KO", "treatment::dmso"],
        )

        expected_counts = self.table.to_dataframe(dense=True).round().astype(int)
        expected_coldata = self.model_metadata.to_dataframe().loc[
            expected_counts.columns, ["genotype", "treatment"]
        ]
        expected_coldata = expected_coldata.astype(str)
        expected_counts.columns.name = observed_counts.columns.name
        expected_coldata.index.name = observed_coldata.index.name

        assert_frame_equal(observed_counts, expected_counts)
        assert_frame_equal(observed_coldata, expected_coldata)
        self.assertEqual(
            observed_formula, "genotype + treatment + genotype:treatment"
        )
        self.assertEqual(observed_reference_levels, ["genotype::KO", "treatment::dmso"])
        self.assertEqual(observed_reduced, "")

    def test_prepare_model_inputs_rejects_missing_formula_columns(self):
        with self.assertRaisesRegex(ValueError, "missing columns required by the model"):
            methods._prepare_model_inputs(
                self.table,
                self.model_metadata,
                fixed_effects_formula="genotype + missing_column",
                min_total_count=0,
            )

    def test_validate_effect_specs_rejects_specs_for_lrt(self):
        with self.assertRaisesRegex(
            ValueError, "effect_specs are not supported when test='lrt'"
        ):
            methods._validate_effect_specs(
                ["contrast::condition::treated::control"], test="lrt"
            )

    def test_prepare_inputs_requires_two_overlapping_samples(self):
        with self.assertRaisesRegex(ValueError, "At least two samples must overlap"):
            methods._prepare_inputs(
                self.table,
                self.condition_unmatched,
                min_total_count=0,
                reference_level="",
            )

    def test_prepare_inputs_requires_two_unique_levels(self):
        with self.assertRaisesRegex(
            ValueError, "Condition metadata must contain at least two unique levels"
        ):
            methods._prepare_inputs(
                self.table,
                self.condition_single_level,
                min_total_count=0,
                reference_level="",
            )

    def test_prepare_inputs_rejects_negative_counts(self):
        negative_counts = self.table.to_dataframe(dense=True)
        negative_counts.iloc[0, 0] = -1
        negative_table = biom.Table(
            negative_counts.values,
            observation_ids=negative_counts.index,
            sample_ids=negative_counts.columns,
        )

        with self.assertRaisesRegex(
            ValueError, "Feature table counts must be non-negative integers"
        ):
            methods._prepare_inputs(
                negative_table,
                self.condition,
                min_total_count=0,
                reference_level="",
            )

    def test_prepare_inputs_rejects_empty_post_filter_table(self):
        with self.assertRaisesRegex(ValueError, "No genes remain after filtering"):
            methods._prepare_inputs(
                self.table,
                self.condition,
                min_total_count=1000,
                reference_level="",
            )

    def test_prepare_inputs_requires_reference_for_more_than_two_levels(self):
        with self.assertRaisesRegex(
            ValueError, "Set reference_level to define the baseline"
        ):
            methods._prepare_inputs(
                self.table,
                self.condition_three_levels,
                min_total_count=0,
                reference_level="",
            )

    def test_prepare_inputs_rejects_missing_reference_level(self):
        with self.assertRaisesRegex(
            ValueError, 'reference_level "missing" is not present'
        ):
            methods._prepare_inputs(
                self.table,
                self.condition,
                min_total_count=0,
                reference_level="missing",
            )

    def test_write_r_script_contains_expected_steps(self):
        with TemporaryDirectory() as temp_dir:
            script_fp = Path(temp_dir) / "run_deseq2.R"

            methods._write_r_script(script_fp)

            script = script_fp.read_text(encoding="utf-8")

        self.assertIn("fixed_effects_formula <- get_arg(\"--fixed-effects-formula\")", script)
        self.assertIn("results_names <- resultsNames(dds)", script)
        self.assertIn("effect_specs <- read_list_file(effect_specs_path)", script)
        self.assertIn("test_kind == \"lrt\"", script)
        self.assertIn("simple_dds_cache <- list()", script)
        self.assertIn("sample_distances_path <- get_arg(\"--sample-distances\")", script)
        self.assertIn(
            "sample_distance_order_path <- get_arg(\"--sample-distance-order\")",
            script,
        )
        self.assertIn("vsd <- vst(dds, blind = FALSE", script)
        self.assertIn("sample_hclust <- hclust(sample_dists)", script)
        self.assertNotIn("--ma-plot", script)
        self.assertNotIn("--volcano-plot", script)
        self.assertNotIn("plotMA(", script)
        self.assertNotIn("--test-level", script)
        self.assertTrue(script.endswith("\n"))

    @patch("q2_deseq2.methods.run")
    def test_run_deseq2_executes_command_and_reads_outputs(self, run_mock):
        expected_results = self.sample_run_result.results
        expected_normalized_counts = self.sample_run_result.normalized_counts
        expected_sample_distances = self.sample_run_result.sample_distance_matrix
        expected_sample_distance_order = self.sample_run_result.sample_distance_order

        def _fake_run(cmd, check, capture_output, text):
            self.assertTrue(check)
            self.assertTrue(capture_output)
            self.assertTrue(text)
            self.assertEqual(cmd[0], "Rscript")

            counts_fp = Path(cmd[cmd.index("--counts") + 1])
            coldata_fp = Path(cmd[cmd.index("--coldata") + 1])
            results_fp = Path(cmd[cmd.index("--results") + 1])
            normalized_counts_fp = Path(cmd[cmd.index("--normalized-counts") + 1])
            sample_distances_fp = Path(cmd[cmd.index("--sample-distances") + 1])
            sample_distance_order_fp = Path(
                cmd[cmd.index("--sample-distance-order") + 1]
            )
            summary_fp = Path(cmd[cmd.index("--summary") + 1])
            results_names_fp = Path(cmd[cmd.index("--results-names") + 1])
            reference_levels_fp = Path(cmd[cmd.index("--reference-levels") + 1])
            effect_specs_fp = Path(cmd[cmd.index("--effect-specs") + 1])

            self.assertTrue(counts_fp.exists())
            self.assertTrue(coldata_fp.exists())
            self.assertEqual(cmd[cmd.index("--fixed-effects-formula") + 1], "condition")
            self.assertEqual(
                reference_levels_fp.read_text(encoding="utf-8"),
                "condition::control\n",
            )
            self.assertIn("contrast::condition::other::control", effect_specs_fp.read_text(encoding="utf-8"))

            expected_results.to_csv(results_fp, sep="\t", index=False)
            expected_normalized_counts.to_csv(
                normalized_counts_fp, sep="\t", index=False
            )
            expected_sample_distances.to_csv(
                sample_distances_fp, sep="\t", index_label="sample_id"
            )
            sample_distance_order_fp.write_text(
                "\n".join(expected_sample_distance_order) + "\n",
                encoding="utf-8",
            )
            results_names_fp.write_text(
                "Intercept\ncondition_treated_vs_control\n",
                encoding="utf-8",
            )
            summary_fp.write_text("summary\n", encoding="utf-8")
            return Mock(returncode=0)

        run_mock.side_effect = _fake_run

        observed = methods.run_deseq2(
            table=self.table,
            condition=self.condition_three_levels,
            reference_level="control",
            min_total_count=3,
            fit_type="mean",
            alpha=0.01,
            cooks_cutoff=False,
            independent_filtering=False,
        )

        run_mock.assert_called_once()
        command = run_mock.call_args.args[0]
        self.assertIn("--fit-type", command)
        self.assertIn("mean", command)
        self.assertIn("--alpha", command)
        self.assertIn("0.01", command)
        self.assertNotIn("--ma-plot", command)
        self.assertNotIn("--volcano-plot", command)
        self.assertIn("false", command)
        assert_frame_equal(observed.results, expected_results)
        assert_frame_equal(observed.normalized_counts, expected_normalized_counts)
        assert_frame_equal(observed.sample_distance_matrix, expected_sample_distances)
        self.assertEqual(observed.sample_distance_order, expected_sample_distance_order)
        self.assertEqual(
            observed.default_effect_id, "contrast::condition::other::control"
        )
        self.assertEqual(observed.available_results_names, ("Intercept", "condition_treated_vs_control"))

    @patch("q2_deseq2.methods.run")
    def test_run_deseq2_model_executes_command_and_reads_outputs(self, run_mock):
        expected_results = pd.DataFrame(
            [
                {
                    "effect_id": "simple::genotype::nonKO::KO|within::treatment::compoundA",
                    "effect_label": "genotype: nonKO vs KO within treatment=compoundA",
                    "effect_kind": "simple",
                    "effect_expression": 'contrast=c("genotype","nonKO","KO"); within="treatment=compoundA"',
                    "comparison": "nonKO vs. KO (treatment=compoundA)",
                    "test_level": "nonKO",
                    "reference_level": "KO",
                    "feature_id": "GG_OTU_1",
                    "baseMean": 25.0,
                    "log2FoldChange": 3.1,
                    "lfcSE": 0.9,
                    "stat": 3.2,
                    "pvalue": 0.02,
                    "padj": 0.04,
                }
            ]
        )
        expected_normalized_counts = self.sample_run_result.normalized_counts
        expected_sample_distances = self.sample_run_result.sample_distance_matrix
        expected_sample_distance_order = self.sample_run_result.sample_distance_order

        def _fake_run(cmd, check, capture_output, text):
            results_fp = Path(cmd[cmd.index("--results") + 1])
            normalized_counts_fp = Path(cmd[cmd.index("--normalized-counts") + 1])
            sample_distances_fp = Path(cmd[cmd.index("--sample-distances") + 1])
            sample_distance_order_fp = Path(
                cmd[cmd.index("--sample-distance-order") + 1]
            )
            summary_fp = Path(cmd[cmd.index("--summary") + 1])
            results_names_fp = Path(cmd[cmd.index("--results-names") + 1])
            reference_levels_fp = Path(cmd[cmd.index("--reference-levels") + 1])
            effect_specs_fp = Path(cmd[cmd.index("--effect-specs") + 1])

            self.assertEqual(
                cmd[cmd.index("--fixed-effects-formula") + 1],
                "genotype + treatment + genotype:treatment",
            )
            self.assertEqual(cmd[cmd.index("--test") + 1], "wald")
            self.assertEqual(
                reference_levels_fp.read_text(encoding="utf-8"),
                "genotype::KO\ntreatment::dmso\n",
            )
            self.assertEqual(
                effect_specs_fp.read_text(encoding="utf-8"),
                "simple::genotype::nonKO::KO|within::treatment::compoundA\n",
            )

            expected_results.to_csv(results_fp, sep="\t", index=False)
            expected_normalized_counts.to_csv(
                normalized_counts_fp, sep="\t", index=False
            )
            expected_sample_distances.to_csv(
                sample_distances_fp, sep="\t", index_label="sample_id"
            )
            sample_distance_order_fp.write_text(
                "\n".join(expected_sample_distance_order) + "\n",
                encoding="utf-8",
            )
            results_names_fp.write_text(
                "Intercept\ngenotype_nonKO_vs_KO\ntreatment_compoundA_vs_dmso\n",
                encoding="utf-8",
            )
            summary_fp.write_text("summary\n", encoding="utf-8")
            return Mock(returncode=0)

        run_mock.side_effect = _fake_run

        observed = methods.run_deseq2_model(
            table=self.table,
            metadata=self.model_metadata,
            fixed_effects_formula="genotype + treatment + genotype:treatment",
            reference_levels=["genotype::KO", "treatment::dmso"],
            effect_specs=[
                "simple::genotype::nonKO::KO|within::treatment::compoundA"
            ],
            min_total_count=0,
        )

        assert_frame_equal(observed.results, expected_results)
        assert_frame_equal(observed.normalized_counts, expected_normalized_counts)
        assert_frame_equal(observed.sample_distance_matrix, expected_sample_distances)
        self.assertEqual(observed.sample_distance_order, expected_sample_distance_order)
        self.assertEqual(
            observed.default_effect_id,
            "simple::genotype::nonKO::KO|within::treatment::compoundA",
        )
        self.assertEqual(
            observed.reference_levels, ("genotype::KO", "treatment::dmso")
        )
        self.assertEqual(
            observed.available_results_names,
            ("Intercept", "genotype_nonKO_vs_KO", "treatment_compoundA_vs_dmso"),
        )
        self.assertNotIn("--ma-plot", run_mock.call_args.args[0])
        self.assertNotIn("--volcano-plot", run_mock.call_args.args[0])

    @patch("q2_deseq2.methods.run")
    def test_run_deseq2_raises_runtime_error_on_failure(self, run_mock):
        run_mock.side_effect = CalledProcessError(
            returncode=1,
            cmd=["Rscript", "run_deseq2.R"],
            stderr="boom",
        )

        with self.assertRaisesRegex(RuntimeError, "DESeq2 command failed"):
            methods.run_deseq2(
                table=self.table,
                condition=self.condition,
                reference_level="control",
            )
