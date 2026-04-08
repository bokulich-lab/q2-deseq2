# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Lesser General Public License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile
from pathlib import Path
from subprocess import CalledProcessError, run

import pandas as pd
from q2_deseq2.utils.frame_utils import (
    _first_value_from_column,
    _unique_non_empty_values,
)
from q2_deseq2.utils.run_data import DESeq2RunResult

from q2_deseq2.utils.analytics import _compute_run_analytics
from q2_deseq2.utils.prep import (
    _normalize_formula,
    _normalize_reference_levels,
    _validate_effect_specs,
)


def _write_lines(path: Path, lines: list[str]) -> None:
    """Write *lines* to *path*, joined by newlines with a trailing newline appended.

    Writes an empty file (zero bytes) when *lines* is empty or all entries
    produce an empty joined string.

    Example::

        >>> from pathlib import Path
        >>> _write_lines(Path("/tmp/refs.txt"), ["treatment::control", "batch::A"])
        # writes "treatment::control\\nbatch::A\\n" to /tmp/refs.txt
    """
    payload = "\n".join(lines)
    if payload:
        payload += "\n"
    path.write_text(payload, encoding="utf-8")


def _write_r_inputs(
    temp_path: Path,
    counts_df: pd.DataFrame,
    coldata_df: pd.DataFrame,
    reference_levels: list[str],
    effect_specs: list[str],
) -> None:
    """Write all DESeq2 input files into *temp_path*.

    Creates four files used by ``run_deseq2.R``:

    * ``counts.tsv`` — feature-by-sample count matrix.
    * ``coldata.tsv`` — per-sample covariate table.
    * ``reference_levels.txt`` — one ``"column::level"`` spec per line.
    * ``effect_specs.txt`` — one effect specification per line.
    """
    counts_df.to_csv(temp_path / "counts.tsv", sep="\t", index_label="feature_id")
    coldata_df.to_csv(temp_path / "coldata.tsv", sep="\t", index_label="sample_id")
    _write_lines(temp_path / "reference_levels.txt", reference_levels)
    _write_lines(temp_path / "effect_specs.txt", effect_specs)


def _build_r_command(
    temp_path: Path,
    script_fp: Path,
    fixed_effects_formula: str,
    test: str,
    reduced_formula: str,
    fit_type: str,
    size_factor_type: str,
    alpha: float,
    cooks_cutoff: bool,
    independent_filtering: bool,
) -> list[str]:
    """Return the ``Rscript`` command list for invoking ``run_deseq2.R``.

    All file paths are derived from *temp_path* (using the names written by
    :func:`_write_r_inputs`).  The returned list is suitable for passing
    directly to :func:`subprocess.run`.
    """
    return [
        "Rscript",
        str(script_fp),
        "--counts",
        str(temp_path / "counts.tsv"),
        "--coldata",
        str(temp_path / "coldata.tsv"),
        "--results",
        str(temp_path / "deseq2_results.tsv"),
        "--size-factors",
        str(temp_path / "size_factors.tsv"),
        "--vst-counts",
        str(temp_path / "vst_counts.tsv"),
        "--summary",
        str(temp_path / "deseq2_summary.txt"),
        "--ma-plot",
        str(temp_path / "ma_plot.png"),
        "--results-names",
        str(temp_path / "results_names.txt"),
        "--reference-levels",
        str(temp_path / "reference_levels.txt"),
        "--effect-specs",
        str(temp_path / "effect_specs.txt"),
        "--fit-type",
        fit_type,
        "--size-factor-type",
        size_factor_type,
        "--alpha",
        str(alpha),
        "--cooks-cutoff",
        str(cooks_cutoff).lower(),
        "--independent-filtering",
        str(independent_filtering).lower(),
        "--fixed-effects-formula",
        fixed_effects_formula,
        "--test",
        test,
        "--reduced-formula",
        reduced_formula,
    ]


def _read_r_outputs(
    temp_path: Path,
    counts_df: pd.DataFrame,
    coldata_df: pd.DataFrame,
    fixed_effects_formula: str,
    reference_levels: list[str],
    test: str,
    reduced_formula: str,
) -> DESeq2RunResult:
    """Parse all DESeq2 output files and return a :class:`DESeq2RunResult`.

    Reads ``deseq2_results.tsv``, ``size_factors.tsv``, ``vst_counts.tsv``,
    and ``results_names.txt`` from *temp_path*, validates the size-factor
    schema, delegates to :func:`_compute_run_analytics`, and assembles the
    result.

    Raises ``RuntimeError`` if ``size_factors.tsv`` is missing required
    columns.
    """
    results_df = pd.read_csv(temp_path / "deseq2_results.tsv", sep="\t")

    size_factors_df = pd.read_csv(temp_path / "size_factors.tsv", sep="\t")
    required_size_factor_columns = {"sample_id", "size_factor"}
    if not required_size_factor_columns <= set(size_factors_df.columns):
        missing_columns = ", ".join(
            sorted(required_size_factor_columns - set(size_factors_df.columns))
        )
        raise RuntimeError(
            "DESeq2 completed but size_factors.tsv is missing required column(s): "
            f"{missing_columns}."
        )
    size_factors = pd.Series(
        size_factors_df["size_factor"].astype(float).to_numpy(),
        index=size_factors_df["sample_id"].astype(str),
    )

    vst_counts = pd.read_csv(temp_path / "vst_counts.tsv", sep="\t", index_col=0)
    vst_counts.index = vst_counts.index.map(str)
    vst_counts.columns = vst_counts.columns.map(str)
    vst_counts.index.name = None
    vst_counts.columns.name = None

    (
        normalized_counts_df,
        sample_distance_matrix,
        sample_distance_order,
        sample_pca_scores,
        sample_pca_percent_variance,
        count_matrix_heatmap,
    ) = _compute_run_analytics(
        counts_df=counts_df,
        size_factors=size_factors,
        vst_counts=vst_counts,
    )

    available_results_names = _unique_non_empty_values(
        (temp_path / "results_names.txt").read_text(encoding="utf-8").splitlines()
    )
    selected_effect_specs = _unique_non_empty_values(
        results_df["effect_id"] if "effect_id" in results_df.columns else ()
    )
    default_effect_id = _first_value_from_column(results_df, "effect_id")

    sample_metadata_df = coldata_df.copy()
    sample_metadata_df.index = sample_metadata_df.index.map(str)
    sample_metadata_df.columns = sample_metadata_df.columns.map(str)

    return DESeq2RunResult(
        results=results_df,
        normalized_counts=normalized_counts_df,
        sample_metadata=sample_metadata_df,
        test_level=_first_value_from_column(results_df, "test_level"),
        reference_level=_first_value_from_column(results_df, "reference_level"),
        default_effect_id=default_effect_id,
        fixed_effects_formula=fixed_effects_formula,
        reference_levels=tuple(reference_levels),
        test=test,
        reduced_formula=reduced_formula,
        available_results_names=available_results_names,
        selected_effect_specs=selected_effect_specs,
        sample_distance_matrix=sample_distance_matrix,
        sample_distance_order=sample_distance_order,
        sample_pca_scores=sample_pca_scores,
        sample_pca_percent_variance=sample_pca_percent_variance,
        count_matrix_heatmap=count_matrix_heatmap,
    )


def _run_deseq2_with_frames(
    counts_df: pd.DataFrame,
    coldata_df: pd.DataFrame,
    fixed_effects_formula: str,
    reference_levels: list[str] | None = None,
    effect_specs: list[str] | None = None,
    test: str = "wald",
    reduced_formula: str = "",
    fit_type: str = "parametric",
    size_factor_type: str = "poscounts",
    alpha: float = 0.05,
    cooks_cutoff: bool = True,
    independent_filtering: bool = True,
) -> DESeq2RunResult:
    """Run the DESeq2 R script and return structured results.

    Writes input DataFrames to a temporary directory, invokes
    ``run_deseq2.R`` via ``Rscript``, reads back all output files, computes
    post-run analytics (normalised counts, PCA, sample distances, heatmap),
    and returns a :class:`~q2_deseq2.utils.run_data.DESeq2RunResult`.

    Raises:
        ValueError: when *test* is ``"lrt"`` and *reduced_formula* is empty,
            or when *effect_specs* are provided for an LRT run.
        RuntimeError: when the R subprocess exits with a non-zero code or
            any expected output file is missing.
    """
    normalized_reference_levels = _normalize_reference_levels(reference_levels)
    normalized_effect_specs = _validate_effect_specs(effect_specs, test)
    normalized_reduced_formula = (
        _normalize_formula(reduced_formula, parameter_name="reduced_formula")
        if str(reduced_formula).strip()
        else ""
    )
    if test == "lrt" and not normalized_reduced_formula:
        raise ValueError("reduced_formula is required when test='lrt'.")

    script_fp = Path(__file__).resolve().parent.parent / "r" / "run_deseq2.R"

    with tempfile.TemporaryDirectory(prefix="q2-deseq2-") as temp_dir:
        temp_path = Path(temp_dir)

        _write_r_inputs(
            temp_path,
            counts_df,
            coldata_df,
            normalized_reference_levels,
            normalized_effect_specs,
        )

        cmd = _build_r_command(
            temp_path,
            script_fp,
            fixed_effects_formula,
            test,
            normalized_reduced_formula,
            fit_type,
            size_factor_type,
            alpha,
            cooks_cutoff,
            independent_filtering,
        )

        try:
            run(cmd, check=True, capture_output=True, text=True)
        except CalledProcessError as exc:
            detail = exc.stderr.strip() or exc.stdout.strip()
            raise RuntimeError(
                f"DESeq2 command failed with exit code {exc.returncode}: {detail}"
            ) from exc

        expected_outputs = {
            "deseq2_results.tsv": temp_path / "deseq2_results.tsv",
            "size_factors.tsv": temp_path / "size_factors.tsv",
            "vst_counts.tsv": temp_path / "vst_counts.tsv",
            "results_names.txt": temp_path / "results_names.txt",
        }
        for expected_name, path in expected_outputs.items():
            if not path.exists():
                raise RuntimeError(
                    f"DESeq2 completed but expected output file "
                    f'"{expected_name}" was not created.'
                )

        return _read_r_outputs(
            temp_path,
            counts_df,
            coldata_df,
            fixed_effects_formula,
            normalized_reference_levels,
            test,
            normalized_reduced_formula,
        )
