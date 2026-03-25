# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import annotations

import json
from pathlib import Path
from typing import NamedTuple

import pandas as pd

from q2_deseq2.types import DESeq2RunDirectoryFormat


class DESeq2RunResult(NamedTuple):
    results: pd.DataFrame
    normalized_counts: pd.DataFrame
    test_level: str = ""
    reference_level: str = ""
    default_effect_id: str = ""
    fixed_effects_formula: str = ""
    reference_levels: tuple[str, ...] = ()
    test: str = "wald"
    reduced_formula: str = ""
    available_results_names: tuple[str, ...] = ()
    selected_effect_specs: tuple[str, ...] = ()
    sample_distance_matrix: pd.DataFrame | None = None
    sample_distance_order: tuple[str, ...] = ()


def _first_non_empty_string(value) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _first_value_from_column(frame: pd.DataFrame, column_name: str) -> str:
    if column_name not in frame.columns:
        return ""
    for value in frame[column_name]:
        normalized = _first_non_empty_string(value)
        if normalized:
            return normalized
    return ""


def _unique_non_empty_values(values) -> tuple[str, ...]:
    ordered_values = []
    seen = set()
    for value in values:
        normalized = _first_non_empty_string(value)
        if normalized and normalized not in seen:
            ordered_values.append(normalized)
            seen.add(normalized)
    return tuple(ordered_values)


def _resolve_default_effect_id(run_result: DESeq2RunResult) -> str:
    if run_result.default_effect_id:
        return run_result.default_effect_id
    if "effect_id" in run_result.results.columns:
        resolved = _first_value_from_column(run_result.results, "effect_id")
        if resolved:
            return resolved
    comparison = _first_value_from_column(run_result.results, "comparison")
    if comparison:
        return f"legacy::{comparison}"
    return ""


def _write_run_result(path: Path, run_result: DESeq2RunResult, alpha: float) -> None:
    path.mkdir(parents=True, exist_ok=True)

    run_result.results.to_csv(path / "deseq2_results.tsv", sep="\t", index=False)
    run_result.normalized_counts.to_csv(
        path / "normalized_counts.tsv", sep="\t", index=False
    )
    if run_result.sample_distance_matrix is not None:
        run_result.sample_distance_matrix.to_csv(
            path / "sample_distances.tsv",
            sep="\t",
            index_label="sample_id",
        )
    if run_result.sample_distance_order:
        (path / "sample_distance_order.txt").write_text(
            "\n".join(run_result.sample_distance_order) + "\n",
            encoding="utf-8",
        )
    (path / "metadata.json").write_text(
        json.dumps(
            {
                "test_level": run_result.test_level,
                "reference_level": run_result.reference_level,
                "default_effect_id": _resolve_default_effect_id(run_result),
                "fixed_effects_formula": run_result.fixed_effects_formula,
                "reference_levels": list(run_result.reference_levels),
                "test": run_result.test,
                "reduced_formula": run_result.reduced_formula,
                "available_results_names": list(run_result.available_results_names),
                "selected_effect_specs": list(
                    run_result.selected_effect_specs
                    or _unique_non_empty_values(
                        run_result.results["effect_id"]
                        if "effect_id" in run_result.results.columns
                        else ()
                    )
                ),
                "alpha": alpha,
            },
            sort_keys=True,
        ),
        encoding="utf-8",
    )


def write_run_result_artifact(
    run_result: DESeq2RunResult, alpha: float
) -> DESeq2RunDirectoryFormat:
    run_data = DESeq2RunDirectoryFormat()
    _write_run_result(Path(str(run_data.path)), run_result=run_result, alpha=alpha)
    return run_data


def _parse_run_results(
    run_data: DESeq2RunDirectoryFormat,
) -> tuple[DESeq2RunResult, float]:
    run_data_path = Path(str(run_data.path))
    metadata = json.loads((run_data_path / "metadata.json").read_text(encoding="utf-8"))
    alpha = float(metadata["alpha"])
    results = pd.read_csv(run_data_path / "deseq2_results.tsv", sep="\t")

    default_effect_id = _first_non_empty_string(metadata.get("default_effect_id"))
    if not default_effect_id and "effect_id" in results.columns:
        default_effect_id = _first_value_from_column(results, "effect_id")

    test_level = _first_non_empty_string(metadata.get("test_level"))
    if not test_level:
        test_level = _first_value_from_column(results, "test_level")

    reference_level = _first_non_empty_string(metadata.get("reference_level"))
    if not reference_level:
        reference_level = _first_value_from_column(results, "reference_level")

    selected_effect_specs = metadata.get("selected_effect_specs") or ()
    if not selected_effect_specs and "effect_id" in results.columns:
        selected_effect_specs = _unique_non_empty_values(results["effect_id"])

    sample_distances_path = run_data_path / "sample_distances.tsv"
    sample_distance_matrix = None
    if sample_distances_path.exists():
        sample_distance_matrix = pd.read_csv(
            sample_distances_path, sep="\t", index_col=0
        )
        sample_distance_matrix.index = sample_distance_matrix.index.map(str)
        sample_distance_matrix.columns = sample_distance_matrix.columns.map(str)
        sample_distance_matrix.index.name = None
        sample_distance_matrix.columns.name = None

    sample_distance_order_path = run_data_path / "sample_distance_order.txt"
    sample_distance_order = ()
    if sample_distance_order_path.exists():
        sample_distance_order = tuple(
            _unique_non_empty_values(
                sample_distance_order_path.read_text(encoding="utf-8").splitlines()
            )
        )
    elif sample_distance_matrix is not None:
        sample_distance_order = tuple(sample_distance_matrix.index.tolist())

    run_result = DESeq2RunResult(
        results=results,
        normalized_counts=pd.read_csv(
            run_data_path / "normalized_counts.tsv", sep="\t"
        ),
        test_level=test_level,
        reference_level=reference_level,
        default_effect_id=default_effect_id,
        fixed_effects_formula=_first_non_empty_string(
            metadata.get("fixed_effects_formula")
        ),
        reference_levels=tuple(
            _unique_non_empty_values(metadata.get("reference_levels") or ())
        ),
        test=_first_non_empty_string(metadata.get("test")) or "wald",
        reduced_formula=_first_non_empty_string(metadata.get("reduced_formula")),
        available_results_names=tuple(
            _unique_non_empty_values(metadata.get("available_results_names") or ())
        ),
        selected_effect_specs=tuple(_unique_non_empty_values(selected_effect_specs)),
        sample_distance_matrix=sample_distance_matrix,
        sample_distance_order=sample_distance_order,
    )

    return run_result, alpha
