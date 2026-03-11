# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import html
from pathlib import Path
import textwrap

import pandas as pd

import biom

from q2_deseq2._deseq2 import run_deseq2


def _render_index_html(output_dir: Path, contrast_label: str, alpha: float) -> None:
    results_fp = output_dir / 'deseq2_results.tsv'
    result_table = pd.read_csv(results_fp, sep='\t')
    preview_columns = [
        column for column in ['feature_id', 'log2FoldChange', 'pvalue', 'padj']
        if column in result_table.columns
    ]
    if preview_columns:
        preview_table = result_table.loc[:, preview_columns].head(25)
    else:
        preview_table = result_table.head(25)

    preview_html = preview_table.to_html(index=False, border=0, classes='preview')

    index_html = textwrap.dedent(
        f"""\
        <!doctype html>
        <html lang="en">
        <head>
          <meta charset="utf-8">
          <title>DESeq2 Differential Expression</title>
          <style>
            body {{ font-family: Arial, sans-serif; margin: 2rem; color: #1f2937; }}
            h1, h2 {{ margin-bottom: 0.5rem; }}
            ul {{ margin-top: 0.25rem; }}
            table {{ border-collapse: collapse; width: 100%; }}
            th, td {{ border: 1px solid #d1d5db; padding: 0.4rem; text-align: left; }}
            th {{ background: #f3f4f6; }}
            code {{ background: #f3f4f6; padding: 0.1rem 0.3rem; }}
          </style>
        </head>
        <body>
          <h1>DESeq2 Differential Expression Report</h1>
          <p><strong>Contrast:</strong> {html.escape(contrast_label)}</p>
          <p><strong>alpha:</strong> {alpha}</p>
          <h2>Output Files</h2>
          <ul>
            <li><a href="deseq2_results.tsv">deseq2_results.tsv</a></li>
            <li><a href="normalized_counts.tsv">normalized_counts.tsv</a></li>
            <li><a href="deseq2_summary.txt">deseq2_summary.txt</a></li>
            <li><a href="ma_plot.png">ma_plot.png</a></li>
            <li><a href="volcano_plot.png">volcano_plot.png</a></li>
            <li><a href="deseq2_stdout.txt">deseq2_stdout.txt</a></li>
            <li><a href="deseq2_stderr.txt">deseq2_stderr.txt</a></li>
          </ul>
          <h2>Top Results Preview</h2>
          {preview_html}
        </body>
        </html>
        """
    )
    (output_dir / 'index.html').write_text(index_html, encoding='utf-8')


def differential_expression(
    output_dir: str,
    table: biom.Table,
    condition: str,
    test_level: str = '',
    reference_level: str = '',
    min_total_count: int = 10,
    fit_type: str = 'parametric',
    alpha: float = 0.05,
    cooks_cutoff: bool = True,
    independent_filtering: bool = True
):
    run_result = run_deseq2(
        table=table,
        condition=condition,
        test_level=test_level,
        reference_level=reference_level,
        min_total_count=min_total_count,
        fit_type=fit_type,
        alpha=alpha,
        cooks_cutoff=cooks_cutoff,
        independent_filtering=independent_filtering
    )

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    run_result.results.to_csv(output_path / 'deseq2_results.tsv', sep='\t', index=False)
    run_result.normalized_counts.to_csv(output_path / 'normalized_counts.tsv', sep='\t', index=False)
    (output_path / 'deseq2_summary.txt').write_text(run_result.summary, encoding='utf-8')
    (output_path / 'ma_plot.png').write_bytes(run_result.ma_plot_png)
    (output_path / 'volcano_plot.png').write_bytes(run_result.volcano_plot_png)
    (output_path / 'deseq2_stdout.txt').write_text(run_result.stdout, encoding='utf-8')
    (output_path / 'deseq2_stderr.txt').write_text(run_result.stderr, encoding='utf-8')

    contrast_label = f'{run_result.test_level} vs {run_result.reference_level}'
    _render_index_html(output_path, contrast_label=contrast_label, alpha=alpha)
