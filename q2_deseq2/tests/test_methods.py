# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock

import biom
import numpy as np
import pandas as pd

from q2_deseq2._methods import _differential_expression_table
from q2_deseq2._pipelines import estimate
from q2_deseq2._visualizers import _load_loci_annotation_table, _visualize


class _FakeMetadataColumn:
    def __init__(self, values):
        self._series = pd.Series(values)

    def to_series(self):
        return self._series.copy()


class _FakeLociInput:
    def __init__(self, path):
        self.path = Path(path)


class DifferentialExpressionTests(unittest.TestCase):
    def _build_table(self) -> biom.Table:
        data = np.array([[120, 132, 298, 310], [7, 5, 18, 21], [55, 49, 62, 60]])
        return biom.Table(
            data,
            observation_ids=["gene-1", "gene-2", "gene-3"],
            sample_ids=["sample-1", "sample-2", "sample-3", "sample-4"],
        )

    def _build_condition(self) -> _FakeMetadataColumn:
        return _FakeMetadataColumn(
            {
                "sample-1": "control",
                "sample-2": "control",
                "sample-3": "treated",
                "sample-4": "treated",
            }
        )

    def test_load_loci_annotation_table(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            gff_path = Path(temp_dir) / "genome-a.gff"
            gff_path.write_text(
                "##gff-version 3\n"
                "contig1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene-1;gene=dnaA;product=Replication protein\n"
                "contig1\tsrc\tCDS\t1\t100\t.\t+\t0\tID=cds-1;Parent=gene-1;product=Replication protein\n",
                encoding="utf-8",
            )

            observed = _load_loci_annotation_table(_FakeLociInput(temp_dir))

            self.assertIn("feature_id", observed.columns)
            self.assertIn("gene_name", observed.columns)
            self.assertIn("product", observed.columns)

            annotation = observed.set_index("feature_id")
            self.assertEqual(annotation.loc["gene-1", "gene_name"], "dnaA")
            self.assertEqual(annotation.loc["gene-1", "product"], "Replication protein")

    @mock.patch("q2_deseq2._deseq2.subprocess.run")
    def test_differential_expression_outputs(self, run_mock):
        def fake_run(cmd, check, capture_output, text):
            self.assertEqual(cmd[0], "Rscript")
            results_fp = Path(cmd[cmd.index("--results") + 1])
            normalized_counts_fp = Path(cmd[cmd.index("--normalized-counts") + 1])
            summary_fp = Path(cmd[cmd.index("--summary") + 1])
            ma_plot_fp = Path(cmd[cmd.index("--ma-plot") + 1])
            volcano_plot_fp = Path(cmd[cmd.index("--volcano-plot") + 1])

            results_fp.write_text(
                "feature_id\tbaseMean\tlog2FoldChange\tpvalue\tpadj\n"
                "gene-1\t200\t1.75\t0.0002\t0.0015\n",
                encoding="utf-8",
            )
            normalized_counts_fp.write_text(
                "feature_id\tsample-1\tsample-2\tsample-3\tsample-4\n"
                "gene-1\t121.1\t130.5\t290.2\t305.6\n",
                encoding="utf-8",
            )
            summary_fp.write_text("summary line\n", encoding="utf-8")
            ma_plot_fp.write_bytes(b"PNG")
            volcano_plot_fp.write_bytes(b"PNG")

            return subprocess.CompletedProcess(
                cmd, 0, stdout="DESeq2 run complete", stderr=""
            )

        run_mock.side_effect = fake_run

        with tempfile.TemporaryDirectory() as output_dir:
            _, run_data = _differential_expression_table(
                table=self._build_table(),
                condition=self._build_condition(),
                min_total_count=5,
                fit_type="parametric",
                alpha=0.05,
            )

            _visualize(output_dir=output_dir, differential_expression_run=run_data)

            expected_files = [
                "deseq2_results.tsv",
                "normalized_counts.tsv",
                "deseq2_summary.txt",
                "ma_plot.png",
                "volcano_plot.png",
                "deseq2_stdout.txt",
                "deseq2_stderr.txt",
                "index.html",
            ]
            for filename in expected_files:
                self.assertTrue(Path(output_dir, filename).exists(), filename)

            index_html = Path(output_dir, "index.html").read_text(encoding="utf-8")
            self.assertIn("DESeq2 Differential Expression Report", index_html)
            self.assertIn("treated vs control", index_html)
            stdout_text = Path(output_dir, "deseq2_stdout.txt").read_text(
                encoding="utf-8"
            )
            self.assertIn("DESeq2 run complete", stdout_text)

    @mock.patch("q2_deseq2._deseq2.subprocess.run")
    def test_differential_expression_method_returns_table(self, run_mock):
        def fake_run(cmd, check, capture_output, text):
            results_fp = Path(cmd[cmd.index("--results") + 1])
            normalized_counts_fp = Path(cmd[cmd.index("--normalized-counts") + 1])
            summary_fp = Path(cmd[cmd.index("--summary") + 1])
            ma_plot_fp = Path(cmd[cmd.index("--ma-plot") + 1])
            volcano_plot_fp = Path(cmd[cmd.index("--volcano-plot") + 1])

            results_fp.write_text(
                "feature_id\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n"
                "gene-1\t220\t2.1\t0.2\t10.5\t1e-5\t0.0003\n",
                encoding="utf-8",
            )
            normalized_counts_fp.write_text(
                "feature_id\tsample-1\tsample-2\tsample-3\tsample-4\n"
                "gene-1\t120.0\t130.0\t300.0\t310.0\n",
                encoding="utf-8",
            )
            summary_fp.write_text("summary line\n", encoding="utf-8")
            ma_plot_fp.write_bytes(b"PNG")
            volcano_plot_fp.write_bytes(b"PNG")

            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

        run_mock.side_effect = fake_run

        observed, run_data = _differential_expression_table(
            table=self._build_table(),
            condition=self._build_condition(),
            min_total_count=0,
            fit_type="parametric",
            alpha=0.05,
        )

        self.assertIn("feature_id", observed.columns)
        self.assertIn("log2FoldChange", observed.columns)
        self.assertEqual(observed.iloc[0]["feature_id"], "gene-1")
        self.assertTrue(Path(str(run_data.path), "deseq2_results.tsv").exists())
        self.assertTrue(Path(str(run_data.path), "metadata.json").exists())
        self.assertFalse(Path(str(run_data.path), "deseq2_stdout.txt").exists())
        self.assertFalse(Path(str(run_data.path), "deseq2_stderr.txt").exists())

    @mock.patch("q2_deseq2._deseq2.subprocess.run")
    def test_differential_expression_subprocess_failure(self, run_mock):
        run_mock.side_effect = subprocess.CalledProcessError(
            2, cmd=["Rscript"], stderr="DESeq2 failed in model fitting"
        )

        with self.assertRaisesRegex(
            RuntimeError,
            "DESeq2 command failed with exit code 2: DESeq2 failed in model fitting",
        ):
            _differential_expression_table(
                table=self._build_table(),
                condition=self._build_condition(),
                min_total_count=0,
            )

    @mock.patch("q2_deseq2._deseq2.subprocess.run")
    def test_differential_expression_with_loci_annotations(self, run_mock):
        def fake_run(cmd, check, capture_output, text):
            results_fp = Path(cmd[cmd.index("--results") + 1])
            normalized_counts_fp = Path(cmd[cmd.index("--normalized-counts") + 1])
            summary_fp = Path(cmd[cmd.index("--summary") + 1])
            ma_plot_fp = Path(cmd[cmd.index("--ma-plot") + 1])
            volcano_plot_fp = Path(cmd[cmd.index("--volcano-plot") + 1])

            results_fp.write_text(
                "feature_id\tbaseMean\tlog2FoldChange\tpvalue\tpadj\n"
                "gene-1\t200\t1.75\t0.0002\t0.0015\n",
                encoding="utf-8",
            )
            normalized_counts_fp.write_text(
                "feature_id\tsample-1\tsample-2\tsample-3\tsample-4\n"
                "gene-1\t121.1\t130.5\t290.2\t305.6\n",
                encoding="utf-8",
            )
            summary_fp.write_text("summary line\n", encoding="utf-8")
            ma_plot_fp.write_bytes(b"PNG")
            volcano_plot_fp.write_bytes(b"PNG")

            return subprocess.CompletedProcess(cmd, 0, stdout="ok", stderr="")

        run_mock.side_effect = fake_run

        with tempfile.TemporaryDirectory() as temp_dir:
            annotation_dir = Path(temp_dir) / "loci"
            annotation_dir.mkdir()
            (annotation_dir / "genome-a.gff").write_text(
                "##gff-version 3\n"
                "contig1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene-1;gene=dnaA;product=Replication protein\n",
                encoding="utf-8",
            )

            output_dir = Path(temp_dir) / "viz"
            _, run_data = _differential_expression_table(
                table=self._build_table(),
                condition=self._build_condition(),
                min_total_count=5,
                fit_type="parametric",
                alpha=0.05,
            )

            _visualize(
                output_dir=str(output_dir),
                differential_expression_run=run_data,
                annotations=_FakeLociInput(annotation_dir),
            )

            annotated_results_fp = output_dir / "deseq2_results_annotated.tsv"
            self.assertTrue(annotated_results_fp.exists())
            annotated_results = annotated_results_fp.read_text(encoding="utf-8")
            self.assertIn("gene_name", annotated_results)
            self.assertIn("dnaA", annotated_results)

            index_html = (output_dir / "index.html").read_text(encoding="utf-8")
            self.assertIn("deseq2_results_annotated.tsv", index_html)

    def test_differential_expression_pipeline_returns_table_and_visualization(self):
        calls = []
        table_output = object()
        run_output = object()
        viz_output = object()

        def _table_action(**kwargs):
            calls.append(("table", kwargs))
            return (table_output, run_output)

        def _visualization_action(**kwargs):
            calls.append(("visualizer", kwargs))
            return (viz_output,)

        class _FakeContext:
            def get_action(self, plugin_name, action_name):
                calls.append(("get_action", plugin_name, action_name))
                if action_name == "_differential_expression_table":
                    return _table_action
                if action_name == "_differential_expression":
                    return _visualization_action
                raise KeyError(action_name)

        observed_table, observed_viz = estimate(
            ctx=_FakeContext(),
            table="table-artifact",
            condition="condition-metadata",
            annotations="annotation-artifact",
            test_level="treated",
            reference_level="control",
            min_total_count=5,
            fit_type="parametric",
            alpha=0.05,
            cooks_cutoff=False,
            independent_filtering=False,
        )

        self.assertIs(observed_table, table_output)
        self.assertIs(observed_viz, viz_output)
        self.assertIn(("get_action", "deseq2", "_differential_expression_table"), calls)
        self.assertIn(("get_action", "deseq2", "_differential_expression"), calls)
        self.assertIn(
            (
                "table",
                {
                    "table": "table-artifact",
                    "condition": "condition-metadata",
                    "test_level": "treated",
                    "reference_level": "control",
                    "min_total_count": 5,
                    "fit_type": "parametric",
                    "alpha": 0.05,
                    "cooks_cutoff": False,
                    "independent_filtering": False,
                },
            ),
            calls,
        )
        self.assertIn(
            (
                "visualizer",
                {
                    "differential_expression_run": run_output,
                    "annotations": "annotation-artifact",
                },
            ),
            calls,
        )
