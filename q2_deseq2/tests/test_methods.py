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
            sample_metadata=pd.DataFrame(
                {"condition": ["control", "treated"]},
                index=["Sample1", "Sample2"],
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
            sample_pca_scores=pd.DataFrame(
                {"PC1": [-2.1, 2.1], "PC2": [0.4, -0.4]},
                index=["Sample1", "Sample2"],
            ),
            sample_pca_percent_variance=(68.2, 21.5),
            count_matrix_heatmap=pd.DataFrame(
                {"Sample2": [7.5, 5.1], "Sample1": [6.8, 4.6]},
                index=["GG_OTU_2", "GG_OTU_1"],
            ),
        )

    def _mock_auxiliary_outputs(
        self, counts_df: pd.DataFrame
    ) -> tuple[pd.Series, pd.DataFrame]:
        sample_ids = counts_df.columns.map(str)
        size_factors = pd.Series(
            [float(index + 1) for index in range(len(sample_ids))],
            index=sample_ids,
        )
        vst_counts = counts_df.astype(float).copy()
        for index, sample_id in enumerate(sample_ids, start=1):
            vst_counts[sample_id] = (vst_counts[sample_id] / index) + (index - 1) * 0.5
        return size_factors, vst_counts

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

    def test_prepare_model_inputs_accepts_dataframe_metadata(self):
        metadata_df = self.model_metadata.to_dataframe()

        (
            observed_counts,
            observed_coldata,
            observed_formula,
            observed_reference_levels,
            observed_reduced,
        ) = methods._prepare_model_inputs(
            self.table,
            metadata_df,
            fixed_effects_formula="genotype + treatment + genotype:treatment",
            min_total_count=0,
            reference_levels=["genotype::KO", "treatment::dmso"],
        )

        expected_counts = self.table.to_dataframe(dense=True).round().astype(int)
        expected_coldata = metadata_df.loc[
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

    def test_compute_run_analytics_uses_python_postprocessing(self):
        counts_df = pd.DataFrame(
            {"Sample1": [2, 6], "Sample2": [4, 10]},
            index=["GG_OTU_1", "GG_OTU_2"],
        )
        size_factors = pd.Series({"Sample1": 2.0, "Sample2": 4.0})
        vst_counts = pd.DataFrame(
            {"Sample1": [1.0, 2.0], "Sample2": [4.0, 6.0]},
            index=["GG_OTU_1", "GG_OTU_2"],
        )

        (
            observed_normalized_counts,
            observed_sample_distances,
            observed_sample_distance_order,
            observed_sample_pca_scores,
            observed_sample_pca_percent_variance,
            observed_count_matrix_heatmap,
        ) = methods._compute_run_analytics(
            counts_df=counts_df,
            size_factors=size_factors,
            vst_counts=vst_counts,
        )

        expected_normalized_counts = pd.DataFrame(
            [
                {"feature_id": "GG_OTU_1", "Sample1": 1.0, "Sample2": 1.0},
                {"feature_id": "GG_OTU_2", "Sample1": 3.0, "Sample2": 2.5},
            ]
        )
        expected_sample_distances = pd.DataFrame(
            [[0.0, 5.0], [5.0, 0.0]],
            index=["Sample1", "Sample2"],
            columns=["Sample1", "Sample2"],
        )
        expected_count_matrix_heatmap = pd.DataFrame(
            {"Sample1": [2.0, 1.0], "Sample2": [6.0, 4.0]},
            index=["GG_OTU_2", "GG_OTU_1"],
        )

        assert_frame_equal(observed_normalized_counts, expected_normalized_counts)
        assert_frame_equal(observed_sample_distances, expected_sample_distances)
        self.assertEqual(observed_sample_distance_order, ("Sample1", "Sample2"))
        self.assertAlmostEqual(abs(observed_sample_pca_scores.loc["Sample1", "PC1"]), 2.5)
        self.assertAlmostEqual(abs(observed_sample_pca_scores.loc["Sample2", "PC1"]), 2.5)
        self.assertAlmostEqual(
            observed_sample_pca_scores.loc["Sample1", "PC1"],
            -observed_sample_pca_scores.loc["Sample2", "PC1"],
        )
        self.assertAlmostEqual(observed_sample_pca_scores["PC2"].abs().max(), 0.0)
        self.assertEqual(observed_sample_pca_percent_variance, (100.0, 0.0))
        assert_frame_equal(observed_count_matrix_heatmap, expected_count_matrix_heatmap)

    def test_r_script_contains_expected_steps(self):
        source_path = Path(methods.__file__).parent / "r" / "run_deseq2.R"
        script = source_path.read_text(encoding="utf-8")

        self.assertIn("fixed_effects_formula <- get_arg(\"--fixed-effects-formula\")", script)
        self.assertIn("results_names <- resultsNames(dds)", script)
        self.assertIn("effect_specs <- read_list_file(effect_specs_path)", script)
        self.assertIn("test_kind == \"lrt\"", script)
        self.assertIn("simple_dds_cache <- list()", script)
        self.assertIn("size_factors_path <- get_arg(\"--size-factors\")", script)
        self.assertIn("vst_counts_path <- get_arg(\"--vst-counts\")", script)
        self.assertIn("size_factor_type <- get_arg(\"--size-factor-type\")", script)
        self.assertIn("vsd <- tryCatch(", script)
        self.assertIn("vst(dds, blind = FALSE", script)
        self.assertIn("varianceStabilizingTransformation(dds, blind = FALSE)", script)
        self.assertIn("size_factors_df <- data.frame(", script)
        self.assertIn("vst_counts <- as.data.frame(assay(vsd))", script)
        self.assertIn("sfType = size_factor_type", script)
        self.assertIn("deseq_with_fit_fallback(", script)
        self.assertNotIn("--normalized-counts", script)
        self.assertNotIn("--sample-distances", script)
        self.assertNotIn("--sample-distance-order", script)
        self.assertNotIn("--sample-pca", script)
        self.assertNotIn("--count-matrix-heatmap", script)
        self.assertIn("--ma-plot", script)
        self.assertNotIn("--volcano-plot", script)
        self.assertNotIn("sample_hclust <- hclust(sample_dists)", script)
        self.assertNotIn("sample_pca_df <- data.frame(", script)
        self.assertNotIn("count_matrix_heatmap_df <- as.data.frame(", script)
        self.assertIn("plotMA(", script)
        self.assertNotIn("--test-level", script)
        self.assertTrue(script.endswith("\n"))

    @patch("q2_deseq2._methodlib.runner.run")
    def test_run_deseq2_executes_command_and_reads_outputs(self, run_mock):
        expected_results = self.sample_run_result.results
        captured = {}

        def _fake_run(cmd, check, capture_output, text):
            self.assertTrue(check)
            self.assertTrue(capture_output)
            self.assertTrue(text)
            self.assertEqual(cmd[0], "Rscript")

            counts_fp = Path(cmd[cmd.index("--counts") + 1])
            coldata_fp = Path(cmd[cmd.index("--coldata") + 1])
            results_fp = Path(cmd[cmd.index("--results") + 1])
            size_factors_fp = Path(cmd[cmd.index("--size-factors") + 1])
            vst_counts_fp = Path(cmd[cmd.index("--vst-counts") + 1])
            summary_fp = Path(cmd[cmd.index("--summary") + 1])
            ma_plot_fp = Path(cmd[cmd.index("--ma-plot") + 1])
            results_names_fp = Path(cmd[cmd.index("--results-names") + 1])
            reference_levels_fp = Path(cmd[cmd.index("--reference-levels") + 1])
            effect_specs_fp = Path(cmd[cmd.index("--effect-specs") + 1])

            self.assertTrue(counts_fp.exists())
            self.assertTrue(coldata_fp.exists())
            captured["counts"] = pd.read_csv(counts_fp, sep="\t", index_col=0)
            captured["counts"].index.name = None
            captured["sample_metadata"] = pd.read_csv(coldata_fp, sep="\t", index_col=0)
            captured["sample_metadata"].index.name = None
            self.assertEqual(cmd[cmd.index("--fixed-effects-formula") + 1], "condition")
            self.assertEqual(cmd[cmd.index("--size-factor-type") + 1], "iterate")
            self.assertEqual(
                reference_levels_fp.read_text(encoding="utf-8"),
                "condition::control\n",
            )
            self.assertIn("contrast::condition::other::control", effect_specs_fp.read_text(encoding="utf-8"))

            captured["size_factors"], captured["vst_counts"] = self._mock_auxiliary_outputs(
                captured["counts"]
            )
            expected_results.to_csv(results_fp, sep="\t", index=False)
            captured["size_factors"].rename_axis("sample_id").reset_index(
                name="size_factor"
            ).to_csv(
                size_factors_fp,
                sep="\t",
                index=False,
            )
            captured["vst_counts"].to_csv(
                vst_counts_fp,
                sep="\t",
                index_label="feature_id",
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
            size_factor_type="iterate",
            alpha=0.01,
            cooks_cutoff=False,
            independent_filtering=False,
        )

        run_mock.assert_called_once()
        command = run_mock.call_args.args[0]
        self.assertIn("--fit-type", command)
        self.assertIn("mean", command)
        self.assertIn("--size-factor-type", command)
        self.assertIn("iterate", command)
        self.assertIn("--alpha", command)
        self.assertIn("0.01", command)
        self.assertIn("--size-factors", command)
        self.assertIn("--vst-counts", command)
        self.assertNotIn("--normalized-counts", command)
        self.assertNotIn("--sample-distances", command)
        self.assertNotIn("--sample-pca", command)
        self.assertIn("--ma-plot", command)
        self.assertNotIn("--volcano-plot", command)
        self.assertIn("false", command)
        assert_frame_equal(observed.results, expected_results)
        (
            expected_normalized_counts,
            expected_sample_distances,
            expected_sample_distance_order,
            expected_sample_pca_scores,
            expected_sample_pca_percent_variance,
            expected_count_matrix_heatmap,
        ) = methods._compute_run_analytics(
            counts_df=captured["counts"],
            size_factors=captured["size_factors"],
            vst_counts=captured["vst_counts"],
        )
        assert_frame_equal(observed.normalized_counts, expected_normalized_counts)
        captured["sample_metadata"].index.name = observed.sample_metadata.index.name
        assert_frame_equal(observed.sample_metadata, captured["sample_metadata"])
        assert_frame_equal(observed.sample_distance_matrix, expected_sample_distances)
        self.assertEqual(observed.sample_distance_order, expected_sample_distance_order)
        assert_frame_equal(observed.sample_pca_scores, expected_sample_pca_scores)
        self.assertEqual(
            observed.sample_pca_percent_variance,
            expected_sample_pca_percent_variance,
        )
        assert_frame_equal(observed.count_matrix_heatmap, expected_count_matrix_heatmap)
        self.assertEqual(
            observed.default_effect_id, "contrast::condition::other::control"
        )
        self.assertEqual(observed.available_results_names, ("Intercept", "condition_treated_vs_control"))

    @patch("q2_deseq2.methods.run_deseq2_model")
    def test_run_deseq2_translates_simple_inputs_to_model_runner(
        self, run_deseq2_model_mock
    ):
        run_deseq2_model_mock.return_value = self.sample_run_result._replace(
            test_level="",
            reference_level="",
        )

        observed = methods.run_deseq2(
            table=self.table,
            condition=self.condition_three_levels,
            reference_level="control",
            min_total_count=3,
            fit_type="mean",
            size_factor_type="iterate",
            alpha=0.01,
            cooks_cutoff=False,
            independent_filtering=False,
        )

        run_deseq2_model_mock.assert_called_once()
        kwargs = run_deseq2_model_mock.call_args.kwargs
        expected_metadata = pd.DataFrame(
            {"condition": self.condition_three_levels.to_series().astype(str)}
        )
        expected_metadata = expected_metadata.loc[:, ["condition"]]
        expected_metadata.index.name = kwargs["metadata"].index.name
        assert_frame_equal(kwargs["metadata"], expected_metadata)
        self.assertIs(kwargs["table"], self.table)
        self.assertEqual(kwargs["fixed_effects_formula"], "condition")
        self.assertEqual(kwargs["reference_levels"], ["condition::control"])
        self.assertEqual(
            kwargs["effect_specs"],
            [
                "contrast::condition::other::control",
                "contrast::condition::treated::control",
            ],
        )
        self.assertEqual(kwargs["test"], "wald")
        self.assertEqual(kwargs["reduced_formula"], "")
        self.assertEqual(kwargs["min_total_count"], 3)
        self.assertEqual(kwargs["fit_type"], "mean")
        self.assertEqual(kwargs["size_factor_type"], "iterate")
        self.assertEqual(kwargs["alpha"], 0.01)
        self.assertFalse(kwargs["cooks_cutoff"])
        self.assertFalse(kwargs["independent_filtering"])
        self.assertEqual(observed.test_level, "")
        self.assertEqual(observed.reference_level, "control")

    @patch("q2_deseq2._methodlib.runner.run")
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
        captured = {}

        def _fake_run(cmd, check, capture_output, text):
            counts_fp = Path(cmd[cmd.index("--counts") + 1])
            coldata_fp = Path(cmd[cmd.index("--coldata") + 1])
            results_fp = Path(cmd[cmd.index("--results") + 1])
            size_factors_fp = Path(cmd[cmd.index("--size-factors") + 1])
            vst_counts_fp = Path(cmd[cmd.index("--vst-counts") + 1])
            summary_fp = Path(cmd[cmd.index("--summary") + 1])
            ma_plot_fp = Path(cmd[cmd.index("--ma-plot") + 1])
            results_names_fp = Path(cmd[cmd.index("--results-names") + 1])
            reference_levels_fp = Path(cmd[cmd.index("--reference-levels") + 1])
            effect_specs_fp = Path(cmd[cmd.index("--effect-specs") + 1])
            captured["counts"] = pd.read_csv(counts_fp, sep="\t", index_col=0)
            captured["counts"].index.name = None
            captured["sample_metadata"] = pd.read_csv(coldata_fp, sep="\t", index_col=0)
            captured["sample_metadata"].index.name = None

            self.assertEqual(
                cmd[cmd.index("--fixed-effects-formula") + 1],
                "genotype + treatment + genotype:treatment",
            )
            self.assertEqual(cmd[cmd.index("--test") + 1], "wald")
            self.assertEqual(cmd[cmd.index("--size-factor-type") + 1], "poscounts")
            self.assertEqual(
                reference_levels_fp.read_text(encoding="utf-8"),
                "genotype::KO\ntreatment::dmso\n",
            )
            self.assertEqual(
                effect_specs_fp.read_text(encoding="utf-8"),
                "simple::genotype::nonKO::KO|within::treatment::compoundA\n",
            )

            captured["size_factors"], captured["vst_counts"] = self._mock_auxiliary_outputs(
                captured["counts"]
            )
            expected_results.to_csv(results_fp, sep="\t", index=False)
            captured["size_factors"].rename_axis("sample_id").reset_index(
                name="size_factor"
            ).to_csv(
                size_factors_fp,
                sep="\t",
                index=False,
            )
            captured["vst_counts"].to_csv(
                vst_counts_fp,
                sep="\t",
                index_label="feature_id",
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
        (
            expected_normalized_counts,
            expected_sample_distances,
            expected_sample_distance_order,
            expected_sample_pca_scores,
            expected_sample_pca_percent_variance,
            expected_count_matrix_heatmap,
        ) = methods._compute_run_analytics(
            counts_df=captured["counts"],
            size_factors=captured["size_factors"],
            vst_counts=captured["vst_counts"],
        )
        assert_frame_equal(observed.normalized_counts, expected_normalized_counts)
        captured["sample_metadata"].index.name = observed.sample_metadata.index.name
        assert_frame_equal(observed.sample_metadata, captured["sample_metadata"])
        assert_frame_equal(observed.sample_distance_matrix, expected_sample_distances)
        self.assertEqual(observed.sample_distance_order, expected_sample_distance_order)
        assert_frame_equal(observed.sample_pca_scores, expected_sample_pca_scores)
        self.assertEqual(
            observed.sample_pca_percent_variance,
            expected_sample_pca_percent_variance,
        )
        assert_frame_equal(observed.count_matrix_heatmap, expected_count_matrix_heatmap)
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
        self.assertIn("--size-factors", run_mock.call_args.args[0])
        self.assertIn("--vst-counts", run_mock.call_args.args[0])
        self.assertIn("--ma-plot", run_mock.call_args.args[0])
        self.assertNotIn("--volcano-plot", run_mock.call_args.args[0])

    @patch("q2_deseq2._methodlib.runner.run")
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

    def test_prepare_model_inputs_drops_samples_with_nan_metadata(self):
        metadata_df = self.model_metadata.to_dataframe()
        metadata_df.loc["Sample3", "treatment"] = None

        observed_counts, observed_coldata, _, _, _ = methods._prepare_model_inputs(
            self.table,
            metadata_df,
            fixed_effects_formula="genotype + treatment",
            min_total_count=0,
        )

        self.assertNotIn("Sample3", observed_coldata.index)
        self.assertNotIn("Sample3", observed_counts.columns)
        self.assertIn("Sample1", observed_coldata.index)
        self.assertEqual(observed_coldata.shape[0], 5)

    def test_prepare_model_inputs_rejects_numeric_reference_level_column(self):
        metadata_df = pd.DataFrame(
            {"genotype": ["KO", "KO", "nonKO", "nonKO", "KO", "nonKO"],
             "depth": [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]},
            index=[f"Sample{i}" for i in range(1, 7)],
        )

        with self.assertRaisesRegex(ValueError, "numeric metadata column"):
            methods._prepare_model_inputs(
                self.table,
                metadata_df,
                fixed_effects_formula="genotype + depth",
                min_total_count=0,
                reference_levels=["depth::10.0"],
            )

    def test_prepare_model_inputs_rejects_reference_level_not_in_observed_data(self):
        with self.assertRaisesRegex(
            ValueError, 'level "missing".*not present in metadata column'
        ):
            methods._prepare_model_inputs(
                self.table,
                self.model_metadata,
                fixed_effects_formula="genotype + treatment",
                min_total_count=0,
                reference_levels=["genotype::missing"],
            )

    def test_normalize_formula_rejects_special_characters(self):
        with self.assertRaisesRegex(ValueError, "may only contain"):
            methods._normalize_formula("condition; DROP TABLE", "formula")

        with self.assertRaisesRegex(ValueError, "may only contain"):
            methods._normalize_formula("condition & batch", "formula")

    def test_normalize_size_factor_type_rejects_invalid_values(self):
        with self.assertRaisesRegex(ValueError, "size_factor_type must be one of"):
            methods._normalize_size_factor_type("median")

        with self.assertRaisesRegex(ValueError, "size_factor_type is required"):
            methods._normalize_size_factor_type("")

    @patch("q2_deseq2.methods.run_deseq2_model")
    def test_run_deseq2_two_levels_auto_infers_reference_and_sets_test_level(
        self, run_deseq2_model_mock
    ):
        run_deseq2_model_mock.return_value = self.sample_run_result._replace(
            test_level="",
            reference_level="",
        )

        observed = methods.run_deseq2(
            table=self.table,
            condition=self.condition,
        )

        kwargs = run_deseq2_model_mock.call_args.kwargs
        self.assertEqual(kwargs["reference_levels"], ["condition::control"])
        self.assertEqual(
            kwargs["effect_specs"],
            ["contrast::condition::treated::control"],
        )
        self.assertEqual(observed.test_level, "treated")
        self.assertEqual(observed.reference_level, "control")

    @patch("q2_deseq2._methodlib.runner.run")
    def test_run_deseq2_model_lrt_passes_correct_arguments(self, run_mock):
        expected_results = pd.DataFrame(
            [
                {
                    "effect_id": "lrt::genotype + treatment::vs::1",
                    "effect_label": "LRT: genotype + treatment vs 1",
                    "effect_kind": "lrt",
                    "effect_expression": "full=~ genotype + treatment; reduced=~ 1",
                    "feature_id": "GG_OTU_1",
                    "baseMean": 25.0,
                    "log2FoldChange": 3.1,
                    "lfcSE": 0.9,
                    "stat": 9.2,
                    "pvalue": 0.01,
                    "padj": 0.04,
                }
            ]
        )
        captured = {}

        def _fake_run(cmd, check, capture_output, text):
            captured["cmd"] = cmd
            counts_fp = Path(cmd[cmd.index("--counts") + 1])
            results_fp = Path(cmd[cmd.index("--results") + 1])
            size_factors_fp = Path(cmd[cmd.index("--size-factors") + 1])
            vst_counts_fp = Path(cmd[cmd.index("--vst-counts") + 1])
            summary_fp = Path(cmd[cmd.index("--summary") + 1])
            results_names_fp = Path(cmd[cmd.index("--results-names") + 1])

            counts_df = pd.read_csv(counts_fp, sep="\t", index_col=0)
            captured["size_factors"], captured["vst_counts"] = (
                self._mock_auxiliary_outputs(counts_df)
            )
            expected_results.to_csv(results_fp, sep="\t", index=False)
            captured["size_factors"].rename_axis("sample_id").reset_index(
                name="size_factor"
            ).to_csv(size_factors_fp, sep="\t", index=False)
            captured["vst_counts"].to_csv(
                vst_counts_fp, sep="\t", index_label="feature_id"
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
            fixed_effects_formula="genotype + treatment",
            test="lrt",
            reduced_formula="genotype",
            min_total_count=0,
        )

        cmd = captured["cmd"]
        self.assertEqual(cmd[cmd.index("--test") + 1], "lrt")
        self.assertEqual(cmd[cmd.index("--reduced-formula") + 1], "genotype")
        self.assertEqual(cmd[cmd.index("--fixed-effects-formula") + 1], "genotype + treatment")
        assert_frame_equal(observed.results, expected_results)
        self.assertEqual(observed.test, "lrt")
        self.assertEqual(observed.reduced_formula, "genotype")

    def test_run_deseq2_with_frames_rejects_lrt_without_reduced_formula(self):
        counts_df = pd.DataFrame(
            {"Sample1": [10, 20], "Sample2": [30, 40]},
            index=["GG_OTU_1", "GG_OTU_2"],
        )
        coldata_df = pd.DataFrame(
            {"condition": ["control", "treated"]},
            index=["Sample1", "Sample2"],
        )

        with self.assertRaisesRegex(
            ValueError, "reduced_formula is required when test='lrt'"
        ):
            methods._run_deseq2_with_frames(
                counts_df=counts_df,
                coldata_df=coldata_df,
                fixed_effects_formula="condition",
                test="lrt",
                reduced_formula="",
            )
