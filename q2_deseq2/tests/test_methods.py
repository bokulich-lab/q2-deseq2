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
        self.sample_run_result = DESeq2RunResult(
            results=pd.DataFrame(
                [
                    {
                        "feature_id": "GG_OTU_2",
                        "baseMean": 120.0,
                        "log2FoldChange": 2.4,
                        "lfcSE": 0.4,
                        "stat": 6.0,
                        "pvalue": 0.002,
                        "padj": 0.01,
                    }
                ]
            ),
            normalized_counts=pd.DataFrame(
                [{"feature_id": "GG_OTU_2", "Sample1": 3.0, "Sample2": 4.0}]
            ),
            ma_plot_png=b"ma-bytes",
            volcano_plot_png=b"volcano-bytes",
            test_level="treated",
            reference_level="control",
        )

    def test_prepare_inputs_infers_two_level_contrast_and_filters_features(self):
        observed_counts, observed_coldata, observed_test, observed_reference = (
            methods._prepare_inputs(
                self.table,
                self.condition,
                min_total_count=3,
                test_level="",
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
        self.assertEqual(observed_test, "treated")
        self.assertEqual(observed_reference, "control")
        self.assertTrue(
            all(dtype.kind in {"i", "u"} for dtype in observed_counts.dtypes)
        )

    def test_prepare_inputs_accepts_explicit_contrast_for_three_level_metadata(self):
        observed_counts, observed_coldata, observed_test, observed_reference = (
            methods._prepare_inputs(
                self.table,
                self.condition_three_levels,
                min_total_count=0,
                test_level="treated",
                reference_level="other",
            )
        )

        expected_counts = self.table.to_dataframe(dense=True).round().astype(int)
        expected_coldata = pd.DataFrame(
            {"condition": self.condition_three_levels.to_series().astype(str)}
        )

        assert_frame_equal(observed_counts, expected_counts)
        assert_frame_equal(observed_coldata, expected_coldata)
        self.assertEqual(observed_test, "treated")
        self.assertEqual(observed_reference, "other")

    def test_prepare_inputs_requires_two_overlapping_samples(self):
        with self.assertRaisesRegex(ValueError, "At least two samples must overlap"):
            methods._prepare_inputs(
                self.table,
                self.condition_unmatched,
                min_total_count=0,
                test_level="",
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
                test_level="",
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
                test_level="",
                reference_level="",
            )

    def test_prepare_inputs_rejects_empty_post_filter_table(self):
        with self.assertRaisesRegex(ValueError, "No genes remain after filtering"):
            methods._prepare_inputs(
                self.table,
                self.condition,
                min_total_count=1000,
                test_level="",
                reference_level="",
            )

    def test_prepare_inputs_requires_both_levels_or_neither(self):
        with self.assertRaisesRegex(
            ValueError, "Provide both test_level and reference_level"
        ):
            methods._prepare_inputs(
                self.table,
                self.condition,
                min_total_count=0,
                test_level="treated",
                reference_level="",
            )

    def test_prepare_inputs_requires_explicit_contrast_for_more_than_two_levels(self):
        with self.assertRaisesRegex(
            ValueError, "Condition metadata has more than two levels"
        ):
            methods._prepare_inputs(
                self.table,
                self.condition_three_levels,
                min_total_count=0,
                test_level="",
                reference_level="",
            )

    def test_prepare_inputs_rejects_missing_requested_level(self):
        with self.assertRaisesRegex(ValueError, 'test_level "missing" is not present'):
            methods._prepare_inputs(
                self.table,
                self.condition,
                min_total_count=0,
                test_level="missing",
                reference_level="control",
            )

    def test_prepare_inputs_rejects_duplicate_contrast_levels(self):
        with self.assertRaisesRegex(
            ValueError, "test_level and reference_level must be different"
        ):
            methods._prepare_inputs(
                self.table,
                self.condition,
                min_total_count=0,
                test_level="treated",
                reference_level="treated",
            )

    def test_write_r_script_contains_expected_steps(self):
        with TemporaryDirectory() as temp_dir:
            script_fp = Path(temp_dir) / "run_deseq2.R"

            methods._write_r_script(script_fp)

            script = script_fp.read_text(encoding="utf-8")

        self.assertIn("DESeqDataSetFromMatrix", script)
        self.assertIn("--normalized-counts", script)
        self.assertIn("plotMA(res, alpha = alpha", script)
        self.assertIn('ylab = "-log10 adjusted p-value"', script)
        self.assertTrue(script.endswith("\n"))

    @patch("q2_deseq2.methods.run")
    def test_run_deseq2_executes_command_and_reads_outputs(self, run_mock):
        expected_results = self.sample_run_result.results
        expected_normalized_counts = self.sample_run_result.normalized_counts
        expected_ma_plot = self.sample_run_result.ma_plot_png
        expected_volcano_plot = self.sample_run_result.volcano_plot_png

        def _fake_run(cmd, check, capture_output, text):
            self.assertTrue(check)
            self.assertTrue(capture_output)
            self.assertTrue(text)
            self.assertEqual(cmd[0], "Rscript")

            counts_fp = Path(cmd[cmd.index("--counts") + 1])
            coldata_fp = Path(cmd[cmd.index("--coldata") + 1])
            results_fp = Path(cmd[cmd.index("--results") + 1])
            normalized_counts_fp = Path(cmd[cmd.index("--normalized-counts") + 1])
            summary_fp = Path(cmd[cmd.index("--summary") + 1])
            ma_plot_fp = Path(cmd[cmd.index("--ma-plot") + 1])
            volcano_plot_fp = Path(cmd[cmd.index("--volcano-plot") + 1])

            self.assertTrue(counts_fp.exists())
            self.assertTrue(coldata_fp.exists())
            self.assertIn("--test-level", cmd)
            self.assertIn("--reference-level", cmd)

            expected_results.to_csv(results_fp, sep="\t", index=False)
            expected_normalized_counts.to_csv(
                normalized_counts_fp, sep="\t", index=False
            )
            summary_fp.write_text("summary\n", encoding="utf-8")
            ma_plot_fp.write_bytes(expected_ma_plot)
            volcano_plot_fp.write_bytes(expected_volcano_plot)
            return Mock(returncode=0)

        run_mock.side_effect = _fake_run

        observed = methods.run_deseq2(
            table=self.table,
            condition=self.condition,
            test_level="treated",
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
        self.assertIn("false", command)
        assert_frame_equal(observed.results, expected_results)
        assert_frame_equal(observed.normalized_counts, expected_normalized_counts)
        self.assertEqual(observed.ma_plot_png, expected_ma_plot)
        self.assertEqual(observed.volcano_plot_png, expected_volcano_plot)
        self.assertEqual(observed.test_level, "treated")
        self.assertEqual(observed.reference_level, "control")

    @patch("q2_deseq2.methods.run")
    def test_run_deseq2_wraps_called_process_error(self, run_mock):
        run_mock.side_effect = CalledProcessError(
            returncode=23, cmd=["Rscript"], stderr="boom"
        )

        with self.assertRaisesRegex(
            RuntimeError, "DESeq2 command failed with exit code 23: boom"
        ):
            methods.run_deseq2(
                table=self.table,
                condition=self.condition,
                test_level="treated",
                reference_level="control",
            )

    @patch("q2_deseq2.methods.run")
    def test_run_deseq2_raises_when_expected_output_is_missing(self, run_mock):
        def _fake_run(cmd, check, capture_output, text):
            results_fp = Path(cmd[cmd.index("--results") + 1])
            normalized_counts_fp = Path(cmd[cmd.index("--normalized-counts") + 1])
            ma_plot_fp = Path(cmd[cmd.index("--ma-plot") + 1])

            self.sample_run_result.results.to_csv(results_fp, sep="\t", index=False)
            self.sample_run_result.normalized_counts.to_csv(
                normalized_counts_fp, sep="\t", index=False
            )
            ma_plot_fp.write_bytes(self.sample_run_result.ma_plot_png)
            return Mock(returncode=0)

        run_mock.side_effect = _fake_run

        with self.assertRaisesRegex(
            RuntimeError,
            'expected output file "volcano_plot.png" was not created',
        ):
            methods.run_deseq2(
                table=self.table,
                condition=self.condition,
                test_level="treated",
                reference_level="control",
            )

    @patch("q2_deseq2.methods.write_run_result_artifact")
    @patch("q2_deseq2.methods.run_deseq2")
    def test_estimate_differential_expression_invokes_run_and_writer(
        self, run_deseq2_mock, write_run_result_artifact_mock
    ):
        run_deseq2_mock.return_value = self.sample_run_result
        write_run_result_artifact_mock.return_value = "artifact"

        observed_results, observed_run_data = methods._estimate_differential_expression(
            table=self.table,
            condition=self.condition,
            test_level="treated",
            reference_level="control",
            min_total_count=11,
            fit_type="local",
            alpha=0.01,
            cooks_cutoff=False,
            independent_filtering=False,
        )

        run_deseq2_mock.assert_called_once_with(
            table=self.table,
            condition=self.condition,
            test_level="treated",
            reference_level="control",
            min_total_count=11,
            fit_type="local",
            alpha=0.01,
            cooks_cutoff=False,
            independent_filtering=False,
        )
        write_run_result_artifact_mock.assert_called_once_with(
            run_result=self.sample_run_result, alpha=0.01
        )
        assert_frame_equal(observed_results, self.sample_run_result.results)
        self.assertEqual(observed_run_data, "artifact")
