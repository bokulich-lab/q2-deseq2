# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json
from pathlib import Path
from typing import NamedTuple

import pandas as pd

from q2_deseq2.types import DESeq2RunDirectoryFormat


class DESeq2RunResult(NamedTuple):
    results: pd.DataFrame
    normalized_counts: pd.DataFrame
    ma_plot_png: bytes
    volcano_plot_png: bytes
    test_level: str
    reference_level: str


def _write_run_result(path: Path, run_result: DESeq2RunResult, alpha: float) -> None:
    path.mkdir(parents=True, exist_ok=True)

    run_result.results.to_csv(path / "deseq2_results.tsv", sep="\t", index=False)
    run_result.normalized_counts.to_csv(
        path / "normalized_counts.tsv", sep="\t", index=False
    )
    (path / "ma_plot.png").write_bytes(run_result.ma_plot_png)
    (path / "volcano_plot.png").write_bytes(run_result.volcano_plot_png)
    (path / "metadata.json").write_text(
        json.dumps(
            {
                "test_level": run_result.test_level,
                "reference_level": run_result.reference_level,
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

    run_result = DESeq2RunResult(
        results=pd.read_csv(run_data_path / "deseq2_results.tsv", sep="\t"),
        normalized_counts=pd.read_csv(
            run_data_path / "normalized_counts.tsv", sep="\t"
        ),
        ma_plot_png=(run_data_path / "ma_plot.png").read_bytes(),
        volcano_plot_png=(run_data_path / "volcano_plot.png").read_bytes(),
        test_level=str(metadata["test_level"]),
        reference_level=str(metadata["reference_level"]),
    )

    return run_result, alpha
