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

from q2_deseq2._methods import differential_expression_table
from q2_deseq2._visualizers import differential_expression


class _FakeMetadataColumn:
    def __init__(self, values):
        self._series = pd.Series(values)

    def to_series(self):
        return self._series.copy()


class DifferentialExpressionTests(unittest.TestCase):
    def _build_table(self) -> biom.Table:
        data = np.array([
            [120, 132, 298, 310],
            [7, 5, 18, 21],
            [55, 49, 62, 60]
        ])
        return biom.Table(
            data,
            observation_ids=['gene-1', 'gene-2', 'gene-3'],
            sample_ids=['sample-1', 'sample-2', 'sample-3', 'sample-4']
        )

    def _build_condition(self) -> _FakeMetadataColumn:
        return _FakeMetadataColumn({
            'sample-1': 'control',
            'sample-2': 'control',
            'sample-3': 'treated',
            'sample-4': 'treated'
        })

    @mock.patch('q2_deseq2._deseq2.subprocess.run')
    def test_differential_expression_outputs(self, run_mock):
        def fake_run(cmd, check, capture_output, text):
            self.assertEqual(cmd[0], 'Rscript')
            results_fp = Path(cmd[cmd.index('--results') + 1])
            normalized_counts_fp = Path(cmd[cmd.index('--normalized-counts') + 1])
            summary_fp = Path(cmd[cmd.index('--summary') + 1])
            ma_plot_fp = Path(cmd[cmd.index('--ma-plot') + 1])

            results_fp.write_text(
                'feature_id\tbaseMean\tlog2FoldChange\tpvalue\tpadj\n'
                'gene-1\t200\t1.75\t0.0002\t0.0015\n',
                encoding='utf-8'
            )
            normalized_counts_fp.write_text(
                'feature_id\tsample-1\tsample-2\tsample-3\tsample-4\n'
                'gene-1\t121.1\t130.5\t290.2\t305.6\n',
                encoding='utf-8'
            )
            summary_fp.write_text('summary line\n', encoding='utf-8')
            ma_plot_fp.write_bytes(b'PNG')

            return subprocess.CompletedProcess(
                cmd,
                0,
                stdout='DESeq2 run complete',
                stderr=''
            )

        run_mock.side_effect = fake_run

        with tempfile.TemporaryDirectory() as output_dir:
            differential_expression(
                output_dir=output_dir,
                table=self._build_table(),
                condition=self._build_condition(),
                min_total_count=5,
                fit_type='parametric',
                alpha=0.05
            )

            expected_files = [
                'deseq2_results.tsv',
                'normalized_counts.tsv',
                'deseq2_summary.txt',
                'ma_plot.png',
                'deseq2_stdout.txt',
                'deseq2_stderr.txt',
                'index.html'
            ]
            for filename in expected_files:
                self.assertTrue(Path(output_dir, filename).exists(), filename)

            index_html = Path(output_dir, 'index.html').read_text(encoding='utf-8')
            self.assertIn('DESeq2 Differential Expression Report', index_html)
            self.assertIn('treated vs control', index_html)

    @mock.patch('q2_deseq2._deseq2.subprocess.run')
    def test_differential_expression_method_returns_table(self, run_mock):
        def fake_run(cmd, check, capture_output, text):
            results_fp = Path(cmd[cmd.index('--results') + 1])
            normalized_counts_fp = Path(cmd[cmd.index('--normalized-counts') + 1])
            summary_fp = Path(cmd[cmd.index('--summary') + 1])
            ma_plot_fp = Path(cmd[cmd.index('--ma-plot') + 1])

            results_fp.write_text(
                'feature_id\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n'
                'gene-1\t220\t2.1\t0.2\t10.5\t1e-5\t0.0003\n',
                encoding='utf-8'
            )
            normalized_counts_fp.write_text(
                'feature_id\tsample-1\tsample-2\tsample-3\tsample-4\n'
                'gene-1\t120.0\t130.0\t300.0\t310.0\n',
                encoding='utf-8'
            )
            summary_fp.write_text('summary line\n', encoding='utf-8')
            ma_plot_fp.write_bytes(b'PNG')

            return subprocess.CompletedProcess(cmd, 0, stdout='', stderr='')

        run_mock.side_effect = fake_run

        observed = differential_expression_table(
            table=self._build_table(),
            condition=self._build_condition(),
            min_total_count=0,
            fit_type='parametric',
            alpha=0.05
        )

        self.assertIn('feature_id', observed.columns)
        self.assertIn('log2FoldChange', observed.columns)
        self.assertEqual(observed.iloc[0]['feature_id'], 'gene-1')

    @mock.patch('q2_deseq2._deseq2.subprocess.run')
    def test_differential_expression_subprocess_failure(self, run_mock):
        run_mock.side_effect = subprocess.CalledProcessError(
            2,
            cmd=['Rscript'],
            stderr='DESeq2 failed in model fitting'
        )

        with tempfile.TemporaryDirectory() as output_dir:
            with self.assertRaisesRegex(
                RuntimeError,
                'DESeq2 command failed with exit code 2: DESeq2 failed in model fitting'
            ):
                differential_expression(
                    output_dir=output_dir,
                    table=self._build_table(),
                    condition=self._build_condition(),
                    min_total_count=0
                )
