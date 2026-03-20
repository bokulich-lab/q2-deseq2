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
from urllib.parse import unquote

import pandas as pd
from q2_types.genome_data import LociDirectoryFormat

from q2_deseq2._run_data import _parse_run_results
from q2_deseq2.methods import DESeq2RunResult
from q2_deseq2.types import DESeq2RunDirectoryFormat

_ASSETS_DIR = Path(__file__).resolve().parent / "assets"
_HTML_TEMPLATE_PATH = _ASSETS_DIR / "deseq2_report.html"
_VEGA_SPEC_TEMPLATE_PATH = _ASSETS_DIR / "volcano_spec.json"


def _value_or_none(value):
    if pd.isna(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _parse_gff3_attributes(raw_attributes: str) -> dict[str, str]:
    parsed = {}
    for entry in raw_attributes.split(";"):
        entry = entry.strip()
        if not entry:
            continue
        if "=" in entry:
            key, value = entry.split("=", 1)
        else:
            key, value = entry, ""
        parsed[key.strip()] = unquote(value.strip())
    return parsed


def _first_non_empty(attributes: dict[str, str], keys: list[str]) -> str | None:
    for key in keys:
        value = attributes.get(key, "").strip()
        if value:
            return value.split(",")[0]
    return None


def _resolve_reference_gff(annotations: LociDirectoryFormat, reference_id: str) -> Path:
    if not reference_id:
        raise ValueError("reference_id is required when gene annotations are provided.")

    file_map = annotations.file_dict()
    gff_fp = file_map.get(reference_id)
    if gff_fp is None:
        available_refs = ", ".join(sorted(file_map))
        raise ValueError(
            f'Reference "{reference_id}" was not found in annotations artifact. '
            f"Available references: {available_refs}"
        )

    return Path(gff_fp)


def _load_loci(gff_path: Path) -> pd.DataFrame:
    annotation_lookup = {}
    with gff_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            attributes = _parse_gff3_attributes(parts[8])
            feature_ids = []
            for key in ["ID", "locus_tag", "Name", "gene", "Parent"]:
                value = attributes.get(key, "").strip()
                if value:
                    feature_ids.extend(v.strip() for v in value.split(",") if v.strip())

            if not feature_ids:
                continue

            gene_name = _first_non_empty(
                attributes, ["gene", "Name", "gene_name", "locus_tag"]
            )
            product = _first_non_empty(attributes, ["product", "description", "Note"])

            for feature_id in feature_ids:
                existing = annotation_lookup.get(
                    feature_id,
                    {"feature_id": feature_id, "gene_name": None, "product": None},
                )
                if existing["gene_name"] is None and gene_name is not None:
                    existing["gene_name"] = gene_name
                if existing["product"] is None and product is not None:
                    existing["product"] = product
                annotation_lookup[feature_id] = existing

    if not annotation_lookup:
        return pd.DataFrame(columns=["feature_id", "gene_name", "product"])

    annotation_table = pd.DataFrame(annotation_lookup.values())
    annotation_table["feature_id"] = annotation_table["feature_id"].astype(str)
    return annotation_table


def _add_annotations(
    result_table: pd.DataFrame, annotation_table: pd.DataFrame
) -> pd.DataFrame:
    merged = result_table.copy()
    merged["feature_id"] = merged["feature_id"].astype(str)
    if not annotation_table.empty:
        merged = merged.merge(annotation_table, how="left", on="feature_id")
    return merged


def _build_volcano_records(result_table: pd.DataFrame) -> list[dict]:
    records = []
    for _, row in result_table.iterrows():
        feature_id = row.get("feature_id")
        if pd.isna(feature_id):
            feature_id = ""

        gene_name = row.get("gene_name")
        if pd.isna(gene_name):
            gene_name = None

        product = row.get("product")
        if pd.isna(product):
            product = None

        records.append(
            {
                "feature_id": str(feature_id),
                "gene_name": None if gene_name is None else str(gene_name),
                "product": None if product is None else str(product),
                "log2FoldChange": _value_or_none(row.get("log2FoldChange")),
                "pvalue": _value_or_none(row.get("pvalue")),
                "padj": _value_or_none(row.get("padj")),
                "baseMean": _value_or_none(row.get("baseMean")),
            }
        )

    return records


def _render_index_html(
    output_dir: Path,
    result_table: pd.DataFrame,
    contrast_label: str,
    alpha: float,
    include_annotated_results_file: bool,
) -> None:
    volcano_records = _build_volcano_records(result_table)

    volcano_spec = json.loads(_VEGA_SPEC_TEMPLATE_PATH.read_text(encoding="utf-8"))
    volcano_spec["signals"][0]["value"] = alpha
    volcano_spec["data"][0]["values"] = volcano_records
    volcano_spec_json = json.dumps(volcano_spec, separators=(",", ":"))

    preview_columns = ["feature_id"]
    for optional_column in ["gene_name", "product"]:
        if (
            optional_column in result_table.columns
            and result_table[optional_column].notna().any()
        ):
            preview_columns.append(optional_column)
    preview_columns.extend(
        column
        for column in ["log2FoldChange", "pvalue", "padj"]
        if column in result_table.columns
    )
    preview_table = result_table.loc[:, preview_columns].head(25)
    preview_html = preview_table.to_html(index=False, border=0, classes="preview")

    additional_file_list = ""
    if include_annotated_results_file:
        additional_file_list = '<li><a href="deseq2_results_annotated.tsv">deseq2_results_annotated.tsv</a></li>'

    index_html_template = _HTML_TEMPLATE_PATH.read_text(encoding="utf-8")
    index_html = (
        index_html_template.replace("__CONTRAST__", html.escape(contrast_label))
        .replace("__ALPHA__", str(alpha))
        .replace("__PREVIEW_TABLE__", preview_html)
        .replace("__ANNOTATED_RESULTS_FILE__", additional_file_list)
        .replace("__VOLCANO_SPEC__", volcano_spec_json)
    )
    (output_dir / "index.html").write_text(index_html, encoding="utf-8")


def _write_visualization_output(
    output_path: Path,
    run_result: DESeq2RunResult,
    alpha: float,
    display_results: pd.DataFrame | None = None,
    include_annotated_results: bool = False,
) -> None:
    run_result.results.to_csv(output_path / "deseq2_results.tsv", sep="\t", index=False)
    if display_results is None:
        display_results = run_result.results
    elif include_annotated_results:
        display_results.to_csv(
            output_path / "deseq2_results_annotated.tsv", sep="\t", index=False
        )

    run_result.normalized_counts.to_csv(
        output_path / "normalized_counts.tsv", sep="\t", index=False
    )
    (output_path / "ma_plot.png").write_bytes(run_result.ma_plot_png)
    (output_path / "volcano_plot.png").write_bytes(run_result.volcano_plot_png)

    contrast_label = f"{run_result.test_level} vs. {run_result.reference_level}"
    _render_index_html(
        output_path,
        result_table=display_results,
        contrast_label=contrast_label,
        alpha=alpha,
        include_annotated_results_file=include_annotated_results,
    )


def _visualize(
    output_dir: str,
    deseq2_results: DESeq2RunDirectoryFormat,
    gene_annotations: LociDirectoryFormat = None,
    reference_id: str = None,
):
    run_result, alpha = _parse_run_results(deseq2_results)
    annotated_results, include_annotated_results = None, False
    if gene_annotations is not None:
        gff_path = _resolve_reference_gff(gene_annotations, reference_id)
        annotation_table = _load_loci(gff_path)
        annotated_results = _add_annotations(run_result.results, annotation_table)
        include_annotated_results = True

    _write_visualization_output(
        Path(output_dir),
        run_result=run_result,
        alpha=alpha,
        display_results=annotated_results,
        include_annotated_results=include_annotated_results,
    )
