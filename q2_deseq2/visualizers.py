# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
from pathlib import Path
from shutil import copytree
from urllib.parse import unquote

import pandas as pd
from q2_types.genome_data import LociDirectoryFormat

from q2_deseq2._run_data import DESeq2RunResult, _parse_run_results
from q2_deseq2.types import DESeq2RunDirectoryFormat

try:
    import q2templates
except ImportError:  # pragma: no cover - exercised in QIIME 2 environments
    q2templates = None

_ASSETS_DIR = Path(__file__).resolve().parent / "assets" / "deseq2"


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


def _build_plot_records(result_table: pd.DataFrame) -> list[dict]:
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


def _load_vega_spec(spec_name: str) -> dict:
    spec_path = _ASSETS_DIR / "vega" / f"{spec_name}.json"
    return json.loads(spec_path.read_text(encoding="utf-8"))


def _format_comparison_label(test_level: str, reference_level: str) -> str:
    return f"{test_level} vs. {reference_level}"


def _ensure_comparison_columns(
    result_table: pd.DataFrame,
    default_test_level: str,
    reference_level: str,
) -> pd.DataFrame:
    enriched = result_table.copy()
    if "test_level" not in enriched.columns:
        enriched["test_level"] = default_test_level
    else:
        enriched["test_level"] = enriched["test_level"].fillna(default_test_level)

    if "reference_level" not in enriched.columns:
        enriched["reference_level"] = reference_level
    else:
        enriched["reference_level"] = enriched["reference_level"].fillna(
            reference_level
        )

    if "comparison" not in enriched.columns:
        enriched["comparison"] = ""
    comparison_missing = enriched["comparison"].isna() | (
        enriched["comparison"].astype(str).str.strip() == ""
    )
    if comparison_missing.any():
        enriched.loc[comparison_missing, "comparison"] = enriched.loc[
            comparison_missing
        ].apply(
            lambda row: _format_comparison_label(
                str(row["test_level"]), str(row["reference_level"])
            ),
            axis=1,
        )

    return enriched


def _legacy_effect_id(
    comparison: str, test_level: str, reference_level: str
) -> str:
    if test_level and reference_level:
        return f"legacy::{test_level}::{reference_level}"
    if comparison:
        return f"legacy::{comparison}"
    return "legacy::effect"


def _ensure_effect_columns(
    result_table: pd.DataFrame,
    default_effect_id: str,
    default_test_level: str,
    reference_level: str,
) -> tuple[pd.DataFrame, str]:
    enriched = _ensure_comparison_columns(
        result_table,
        default_test_level=default_test_level,
        reference_level=reference_level,
    )

    if "effect_id" not in enriched.columns:
        enriched["effect_id"] = ""
    if "effect_label" not in enriched.columns:
        enriched["effect_label"] = ""
    if "effect_kind" not in enriched.columns:
        enriched["effect_kind"] = ""
    if "effect_expression" not in enriched.columns:
        enriched["effect_expression"] = ""

    for index, row in enriched.iterrows():
        comparison = str(row.get("comparison", "") or "").strip()
        test_level = str(row.get("test_level", "") or "").strip()
        reference = str(row.get("reference_level", "") or "").strip()
        effect_id = str(row.get("effect_id", "") or "").strip()
        effect_label = str(row.get("effect_label", "") or "").strip()
        effect_kind = str(row.get("effect_kind", "") or "").strip()
        effect_expression = str(row.get("effect_expression", "") or "").strip()

        if not effect_id:
            effect_id = _legacy_effect_id(comparison, test_level, reference)
            enriched.at[index, "effect_id"] = effect_id
        if not effect_label:
            effect_label = comparison or effect_id
            enriched.at[index, "effect_label"] = effect_label
        if not effect_kind:
            enriched.at[index, "effect_kind"] = "comparison"
        if not effect_expression:
            enriched.at[index, "effect_expression"] = comparison or effect_label

    available_effect_ids = [
        str(effect_id).strip()
        for effect_id in enriched["effect_id"].dropna().tolist()
        if str(effect_id).strip()
    ]
    if default_effect_id not in available_effect_ids:
        default_effect_id = available_effect_ids[0] if available_effect_ids else ""

    return enriched, default_effect_id


def _collect_effect_options(
    result_table: pd.DataFrame, default_effect_id: str
) -> tuple[list[dict[str, str]], str, str]:
    effect_frame = (
        result_table.loc[:, ["effect_id", "effect_label", "effect_kind"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    effect_options = effect_frame.to_dict(orient="records")
    if not effect_options:
        effect_options = [
            {
                "effect_id": default_effect_id or "legacy::effect",
                "effect_label": "Effect",
                "effect_kind": "comparison",
            }
        ]

    known_effect_ids = {option["effect_id"] for option in effect_options}
    if default_effect_id not in known_effect_ids:
        default_effect_id = effect_options[0]["effect_id"]

    effect_options.sort(key=lambda option: option["effect_id"] != default_effect_id)
    default_effect_label = next(
        option["effect_label"]
        for option in effect_options
        if option["effect_id"] == default_effect_id
    )
    return effect_options, default_effect_id, default_effect_label


def _summarize_results(
    result_table: pd.DataFrame, alpha: float
) -> dict[str, int | str]:
    log2fc = pd.to_numeric(result_table.get("log2FoldChange"), errors="coerce")
    padj = pd.to_numeric(result_table.get("padj"), errors="coerce")
    plottable_mask = log2fc.notna() & padj.notna() & (padj > 0)
    significant_mask = plottable_mask & (padj <= alpha)

    summary = {
        "total_features": int(len(result_table)),
        "plottable_features": int(plottable_mask.sum()),
        "significant_features": int(significant_mask.sum()),
        "alpha_display": str(alpha),
    }
    if {"gene_name", "product"} & set(result_table.columns):
        annotated_mask = pd.Series(False, index=result_table.index)
        for column in ["gene_name", "product"]:
            if column in result_table.columns:
                annotated_mask = annotated_mask | result_table[column].notna()
        summary["annotated_features"] = int(annotated_mask.sum())

    return summary


def _prepare_table_payload(result_table: pd.DataFrame) -> str:
    serializable = result_table.astype(object).where(pd.notna(result_table), None)
    columns = [str(column) for column in serializable.columns]
    rows = serializable.loc[:, columns].values.tolist()
    return json.dumps({"columns": columns, "data": rows}).replace("NaN", "null")


def _copy_report_assets(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for folder in ["css", "js", "vega"]:
        src_dir = _ASSETS_DIR / folder
        if src_dir.exists():
            copytree(src_dir, output_dir / folder, dirs_exist_ok=True)


def _render_report(
    output_dir: Path,
    result_table: pd.DataFrame,
    default_effect_id: str,
    default_test_level: str,
    reference_level: str,
    alpha: float,
    include_annotated_results_file: bool,
) -> None:
    if q2templates is None:
        raise ImportError("q2templates is required to render the DESeq2 visualization.")

    _copy_report_assets(output_dir)
    data_dir = output_dir / "data"
    data_dir.mkdir(exist_ok=True)
    report_results, default_effect_id = _ensure_effect_columns(
        result_table,
        default_effect_id=default_effect_id,
        default_test_level=default_test_level,
        reference_level=reference_level,
    )
    (data_dir / "results_table.json").write_text(
        _prepare_table_payload(report_results), encoding="utf-8"
    )

    volcano_spec = _load_vega_spec("volcano")
    volcano_spec["signals"][0]["value"] = alpha
    ma_spec = _load_vega_spec("ma")
    ma_spec["signals"][0]["value"] = alpha
    effect_options, default_effect_id, default_effect_label = _collect_effect_options(
        report_results, default_effect_id
    )
    default_results = report_results.loc[
        report_results["effect_id"] == default_effect_id
    ].reset_index(drop=True)
    if default_results.empty:
        default_results = report_results.reset_index(drop=True)

    summary = _summarize_results(default_results, alpha)

    templates = [
        str(_ASSETS_DIR / "index.html"),
        str(_ASSETS_DIR / "table.html"),
    ]
    context = {
        "tabs": [
            {"title": "Overview", "url": "index.html"},
            {"title": "Results table", "url": "table.html"},
        ],
        "alpha": str(alpha),
        "summary": summary,
        "vega_volcano_spec": json.dumps(volcano_spec),
        "vega_ma_spec": json.dumps(ma_spec),
        "results_data_path_json": json.dumps("data/results_table.json"),
        "effect_options_json": json.dumps(effect_options),
        "default_effect_id_json": json.dumps(default_effect_id),
        "default_effect_label": default_effect_label,
        "include_annotated_results_file": include_annotated_results_file,
    }
    q2templates.render(templates, str(output_dir), context=context)


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

    _render_report(
        output_path,
        result_table=display_results,
        default_effect_id=run_result.default_effect_id,
        default_test_level=run_result.test_level,
        reference_level=run_result.reference_level,
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
