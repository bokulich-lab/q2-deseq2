# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import html
import json
from pathlib import Path
import textwrap
from urllib.parse import unquote

import biom
import pandas as pd
from q2_types.genome_data import LociDirectoryFormat

from q2_deseq2._deseq2 import DESeq2RunResult, run_deseq2


def _value_or_none(value):
    if pd.isna(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _parse_gff3_attributes(raw_attributes: str) -> dict[str, str]:
    parsed = {}
    for entry in raw_attributes.split(';'):
        entry = entry.strip()
        if not entry:
            continue
        if '=' in entry:
            key, value = entry.split('=', 1)
        else:
            key, value = entry, ''
        parsed[key.strip()] = unquote(value.strip())
    return parsed


def _first_non_empty(attributes: dict[str, str], keys: list[str]) -> str | None:
    for key in keys:
        value = attributes.get(key, '').strip()
        if value:
            return value.split(',')[0]
    return None


def _load_loci_annotation_table(annotations) -> pd.DataFrame:
    loci_path = getattr(annotations, 'path', annotations)
    loci_path = Path(str(loci_path))

    if loci_path.is_file() and loci_path.suffix == '.gff':
        gff_paths = [loci_path]
    else:
        gff_paths = sorted(loci_path.glob('*.gff'))

    if not gff_paths:
        return pd.DataFrame(columns=['feature_id', 'gene_name', 'product'])

    annotation_lookup = {}
    for gff_path in gff_paths:
        with gff_path.open('r', encoding='utf-8') as handle:
            for line in handle:
                if not line.strip() or line.startswith('#'):
                    continue

                parts = line.rstrip('\n').split('\t')
                if len(parts) != 9:
                    continue

                attributes = _parse_gff3_attributes(parts[8])
                feature_ids = []
                for key in ['ID', 'locus_tag', 'Name', 'gene', 'Parent']:
                    value = attributes.get(key, '').strip()
                    if value:
                        feature_ids.extend(v.strip() for v in value.split(',') if v.strip())

                if not feature_ids:
                    continue

                gene_name = _first_non_empty(attributes, ['gene', 'Name', 'gene_name', 'locus_tag'])
                product = _first_non_empty(attributes, ['product', 'description', 'Note'])

                for feature_id in feature_ids:
                    record = annotation_lookup.get(feature_id, {
                        'feature_id': feature_id,
                        'gene_name': None,
                        'product': None
                    })
                    if record['gene_name'] is None and gene_name is not None:
                        record['gene_name'] = gene_name
                    if record['product'] is None and product is not None:
                        record['product'] = product
                    annotation_lookup[feature_id] = record

    if not annotation_lookup:
        return pd.DataFrame(columns=['feature_id', 'gene_name', 'product'])

    annotation_table = pd.DataFrame(annotation_lookup.values())
    annotation_table['feature_id'] = annotation_table['feature_id'].astype(str)
    return annotation_table


def _merge_results_with_annotations(
    result_table: pd.DataFrame, annotation_table: pd.DataFrame
) -> pd.DataFrame:
    merged = result_table.copy()
    merged['feature_id'] = merged['feature_id'].astype(str)
    if not annotation_table.empty:
        merged = merged.merge(annotation_table, how='left', on='feature_id')
    return merged


def _build_volcano_records(result_table: pd.DataFrame) -> list[dict]:
    records = []
    for _, row in result_table.iterrows():
        feature_id = row.get('feature_id')
        if pd.isna(feature_id):
            feature_id = ''

        gene_name = row.get('gene_name')
        if pd.isna(gene_name):
            gene_name = None

        product = row.get('product')
        if pd.isna(product):
            product = None

        records.append({
            'feature_id': str(feature_id),
            'gene_name': None if gene_name is None else str(gene_name),
            'product': None if product is None else str(product),
            'log2FoldChange': _value_or_none(row.get('log2FoldChange')),
            'pvalue': _value_or_none(row.get('pvalue')),
            'padj': _value_or_none(row.get('padj')),
            'baseMean': _value_or_none(row.get('baseMean'))
        })

    return records


def _render_index_html(output_dir: Path, result_table: pd.DataFrame, contrast_label: str,
                       alpha: float, include_annotated_results_file: bool) -> None:
    volcano_records = _build_volcano_records(result_table)
    volcano_spec = {
        '$schema': 'https://vega.github.io/schema/vega/v5.json',
        'description': 'Interactive DESeq2 volcano plot.',
        'width': 860,
        'height': 500,
        'padding': {'left': 75, 'top': 10, 'right': 20, 'bottom': 60},
        'signals': [
            {
                'name': 'alphaCutoff',
                'value': alpha,
                'bind': {
                    'input': 'range',
                    'min': 0.0001,
                    'max': 0.25,
                    'step': 0.0001,
                    'name': 'Adjusted p-value cutoff '
                }
            },
            {
                'name': 'lfcCutoff',
                'value': 1.0,
                'bind': {
                    'input': 'range',
                    'min': 0.0,
                    'max': 5.0,
                    'step': 0.1,
                    'name': '|log2FC| cutoff '
                }
            }
        ],
        'data': [
            {
                'name': 'points',
                'values': volcano_records,
                'transform': [
                    {
                        'type': 'formula',
                        'as': 'plot_p',
                        'expr': 'datum.padj > 0 ? datum.padj : (datum.pvalue > 0 ? datum.pvalue : null)'
                    },
                    {
                        'type': 'formula',
                        'as': 'plot_p_source',
                        'expr': "datum.padj > 0 ? 'padj' : (datum.pvalue > 0 ? 'pvalue' : 'none')"
                    },
                    {
                        'type': 'filter',
                        'expr': ('isValid(datum.log2FoldChange) && isFinite(datum.log2FoldChange) && '
                                 'datum.plot_p != null && isFinite(datum.plot_p)')
                    },
                    {'type': 'formula', 'as': 'neg_log10_p', 'expr': '-log(datum.plot_p)/LN10'},
                    {'type': 'formula', 'as': 'abs_lfc', 'expr': 'abs(datum.log2FoldChange)'},
                    {
                        'type': 'formula',
                        'as': 'sig_flag',
                        'expr': 'datum.plot_p <= alphaCutoff && datum.abs_lfc >= lfcCutoff ? 1 : 0'
                    },
                    {
                        'type': 'formula',
                        'as': 'sig_label',
                        'expr': "datum.sig_flag ? 'Significant' : 'Not significant'"
                    }
                ]
            },
            {
                'name': 'summary',
                'source': 'points',
                'transform': [
                    {
                        'type': 'aggregate',
                        'ops': ['count', 'sum'],
                        'fields': [None, 'sig_flag'],
                        'as': ['total', 'significant']
                    }
                ]
            }
        ],
        'scales': [
            {
                'name': 'x',
                'type': 'linear',
                'domain': {'data': 'points', 'field': 'log2FoldChange'},
                'range': 'width',
                'nice': True,
                'zero': False
            },
            {
                'name': 'y',
                'type': 'linear',
                'domain': {'data': 'points', 'field': 'neg_log10_p'},
                'range': 'height',
                'nice': True,
                'zero': True
            },
            {
                'name': 'color',
                'type': 'ordinal',
                'domain': ['Not significant', 'Significant'],
                'range': ['#9ca3af', '#dc2626']
            }
        ],
        'axes': [
            {'orient': 'bottom', 'scale': 'x', 'title': 'log2 fold change'},
            {'orient': 'left', 'scale': 'y', 'title': '-log10 adjusted p-value (fallback to p-value)'}
        ],
        'legends': [
            {'fill': 'color', 'title': 'Classification', 'orient': 'top-right'}
        ],
        'marks': [
            {
                'type': 'rule',
                'encode': {
                    'enter': {
                        'stroke': {'value': '#dc2626'},
                        'strokeDash': {'value': [5, 4]},
                        'strokeWidth': {'value': 1}
                    },
                    'update': {
                        'x': {'value': 0},
                        'x2': {'signal': 'width'},
                        'y': {'scale': 'y', 'signal': '-log(alphaCutoff)/LN10'}
                    }
                }
            },
            {
                'type': 'rule',
                'encode': {
                    'enter': {
                        'stroke': {'value': '#2563eb'},
                        'strokeDash': {'value': [4, 4]},
                        'strokeWidth': {'value': 1}
                    },
                    'update': {
                        'x': {'scale': 'x', 'signal': 'lfcCutoff'},
                        'x2': {'scale': 'x', 'signal': 'lfcCutoff'},
                        'y': {'value': 0},
                        'y2': {'signal': 'height'}
                    }
                }
            },
            {
                'type': 'rule',
                'encode': {
                    'enter': {
                        'stroke': {'value': '#2563eb'},
                        'strokeDash': {'value': [4, 4]},
                        'strokeWidth': {'value': 1}
                    },
                    'update': {
                        'x': {'scale': 'x', 'signal': '-lfcCutoff'},
                        'x2': {'scale': 'x', 'signal': '-lfcCutoff'},
                        'y': {'value': 0},
                        'y2': {'signal': 'height'}
                    }
                }
            },
            {
                'type': 'symbol',
                'from': {'data': 'points'},
                'encode': {
                    'enter': {
                        'size': {'value': 45},
                        'opacity': {'value': 0.75},
                        'stroke': {'value': '#ffffff'},
                        'strokeWidth': {'value': 0.4}
                    },
                    'update': {
                        'x': {'scale': 'x', 'field': 'log2FoldChange'},
                        'y': {'scale': 'y', 'field': 'neg_log10_p'},
                        'fill': {'scale': 'color', 'field': 'sig_label'},
                        'tooltip': {
                            'signal': (
                                "{'gene_id': datum.feature_id, "
                                "'gene_name': datum.gene_name, "
                                "'product': datum.product, "
                                "'log2FoldChange': format(datum.log2FoldChange, '.4f'), "
                                "'plot_p': datum.plot_p, "
                                "'plot_p_source': datum.plot_p_source, "
                                "'padj': datum.padj, "
                                "'pvalue': datum.pvalue, "
                                "'baseMean': datum.baseMean}"
                            )
                        }
                    },
                    'hover': {
                        'opacity': {'value': 1},
                        'size': {'value': 100}
                    }
                }
            },
            {
                'type': 'text',
                'encode': {
                    'enter': {
                        'x': {'value': 8},
                        'y': {'value': -4},
                        'fontSize': {'value': 12},
                        'fill': {'value': '#111827'}
                    },
                    'update': {
                        'text': {
                            'signal': (
                                "data('summary').length ? "
                                "'Significant genes: ' + format(data('summary')[0].significant, ',') + "
                                "' / ' + format(data('summary')[0].total, ',') : "
                                "'No plottable points in result table.'"
                            )
                        }
                    }
                }
            }
        ]
    }
    volcano_spec_json = json.dumps(volcano_spec, separators=(',', ':'))

    preview_columns = ['feature_id']
    for optional_column in ['gene_name', 'product']:
        if optional_column in result_table.columns and result_table[optional_column].notna().any():
            preview_columns.append(optional_column)
    preview_columns.extend(
        column for column in ['log2FoldChange', 'pvalue', 'padj']
        if column in result_table.columns
    )
    preview_table = result_table.loc[:, preview_columns].head(25)
    preview_html = preview_table.to_html(index=False, border=0, classes='preview')

    additional_file_list = ''
    if include_annotated_results_file:
        additional_file_list = (
            '<li><a href="deseq2_results_annotated.tsv">deseq2_results_annotated.tsv</a></li>'
        )

    index_html_template = textwrap.dedent(
        """\
        <!doctype html>
        <html lang="en">
        <head>
          <meta charset="utf-8">
          <title>DESeq2 Differential Expression</title>
          <script src="https://cdn.jsdelivr.net/npm/vega@5"></script>
          <script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>
          <style>
            body { font-family: "Segoe UI", Tahoma, sans-serif; margin: 1.5rem 2rem; color: #1f2937; background: #f9fafb; }
            h1, h2 { margin-bottom: 0.4rem; }
            .header { background: #ffffff; border: 1px solid #e5e7eb; border-radius: 10px; padding: 1rem 1.2rem; margin-bottom: 1rem; }
            .section { background: #ffffff; border: 1px solid #e5e7eb; border-radius: 10px; padding: 1rem 1.2rem; margin-bottom: 1rem; }
            .meta { display: flex; flex-wrap: wrap; gap: 1.2rem; font-size: 0.95rem; }
            .meta code { background: #eff6ff; color: #1d4ed8; border: 1px solid #bfdbfe; padding: 0.1rem 0.35rem; border-radius: 4px; }
            ul { margin-top: 0.3rem; line-height: 1.6; }
            #volcano-view { min-height: 560px; }
            #volcano-fallback { margin-top: 0.8rem; }
            #volcano-fallback img { max-width: 100%; border: 1px solid #d1d5db; border-radius: 6px; }
            #volcano-error { color: #b91c1c; margin-top: 0.5rem; white-space: pre-wrap; }
            table { border-collapse: collapse; width: 100%; font-size: 0.92rem; }
            th, td { border: 1px solid #d1d5db; padding: 0.4rem; text-align: left; }
            th { background: #eef2ff; color: #1e3a8a; }
            .vega-bindings { display: flex; flex-wrap: wrap; gap: 1.3rem; margin-bottom: 0.7rem; font-size: 0.92rem; }
            .vega-bind label { display: flex; align-items: center; gap: 0.45rem; }
            .vega-bind input[type="range"] { width: 240px; }
          </style>
        </head>
        <body>
          <div class="header">
            <h1>DESeq2 Differential Expression Report</h1>
            <div class="meta">
              <div><strong>Contrast:</strong> <code>__CONTRAST__</code></div>
              <div><strong>Default alpha:</strong> <code>__ALPHA__</code></div>
            </div>
          </div>

          <div class="section">
            <h2>Interactive Volcano Plot</h2>
            <p>Hover points to inspect gene IDs and adjust the displayed significance cutoffs with sliders.</p>
            <div id="volcano-view"></div>
            <div id="volcano-error"></div>
            <div id="volcano-fallback">
              <p>Fallback static volcano plot:</p>
              <img src="volcano_plot.png" alt="Static volcano plot">
            </div>
          </div>

          <div class="section">
            <h2>Output Files</h2>
            <ul>
              <li><a href="deseq2_results.tsv">deseq2_results.tsv</a></li>
              __ANNOTATED_RESULTS_FILE__
              <li><a href="normalized_counts.tsv">normalized_counts.tsv</a></li>
              <li><a href="deseq2_summary.txt">deseq2_summary.txt</a></li>
              <li><a href="ma_plot.png">ma_plot.png</a></li>
              <li><a href="volcano_plot.png">volcano_plot.png</a></li>
              <li><a href="deseq2_stdout.txt">deseq2_stdout.txt</a></li>
              <li><a href="deseq2_stderr.txt">deseq2_stderr.txt</a></li>
            </ul>
          </div>

          <div class="section">
            <h2>Top Results Preview</h2>
            __PREVIEW_TABLE__
          </div>

          <script>
            const volcanoSpec = __VOLCANO_SPEC__;
            const volcanoFallback = document.getElementById('volcano-fallback');
            const volcanoError = document.getElementById('volcano-error');
            vegaEmbed('#volcano-view', volcanoSpec, {actions: false, renderer: 'canvas'})
              .then(function() {
                if (volcanoFallback) {
                  volcanoFallback.style.display = 'none';
                }
              })
              .catch(function(err) {
                if (volcanoError) {
                  volcanoError.textContent =
                    'Interactive volcano plot could not be rendered. ' +
                    'The static fallback image is still available below.\\n\\n' + err;
                }
                if (volcanoFallback) {
                  volcanoFallback.style.display = 'block';
                }
              });
          </script>
        </body>
        </html>
        """
    )
    index_html = (index_html_template
                  .replace('__CONTRAST__', html.escape(contrast_label))
                  .replace('__ALPHA__', str(alpha))
                  .replace('__PREVIEW_TABLE__', preview_html)
                  .replace('__ANNOTATED_RESULTS_FILE__', additional_file_list)
                  .replace('__VOLCANO_SPEC__', volcano_spec_json))
    (output_dir / 'index.html').write_text(index_html, encoding='utf-8')


def _write_visualization_output(
    output_path: Path,
    run_result: DESeq2RunResult,
    alpha: float,
    display_results: pd.DataFrame | None = None,
    include_annotated_results_file: bool = False
) -> None:
    output_path.mkdir(parents=True, exist_ok=True)

    run_result.results.to_csv(output_path / 'deseq2_results.tsv', sep='\t', index=False)
    if display_results is None:
        display_results = run_result.results
    elif include_annotated_results_file:
        display_results.to_csv(output_path / 'deseq2_results_annotated.tsv', sep='\t', index=False)

    run_result.normalized_counts.to_csv(output_path / 'normalized_counts.tsv', sep='\t', index=False)
    (output_path / 'deseq2_summary.txt').write_text(run_result.summary, encoding='utf-8')
    (output_path / 'ma_plot.png').write_bytes(run_result.ma_plot_png)
    (output_path / 'volcano_plot.png').write_bytes(run_result.volcano_plot_png)
    (output_path / 'deseq2_stdout.txt').write_text(run_result.stdout, encoding='utf-8')
    (output_path / 'deseq2_stderr.txt').write_text(run_result.stderr, encoding='utf-8')

    contrast_label = f'{run_result.test_level} vs {run_result.reference_level}'
    _render_index_html(
        output_path,
        result_table=display_results,
        contrast_label=contrast_label,
        alpha=alpha,
        include_annotated_results_file=include_annotated_results_file
    )


def differential_expression(
    output_dir: str,
    table: biom.Table,
    condition: str,
    annotations: LociDirectoryFormat =None,
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
    if annotations is None:
        _write_visualization_output(output_path, run_result=run_result, alpha=alpha)
        return

    annotation_table = _load_loci_annotation_table(annotations)
    annotated_results = _merge_results_with_annotations(run_result.results, annotation_table)
    _write_visualization_output(
        output_path,
        run_result=run_result,
        alpha=alpha,
        display_results=annotated_results,
        include_annotated_results_file=True
    )
