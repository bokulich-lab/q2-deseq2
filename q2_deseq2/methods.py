# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from subprocess import run, CalledProcessError
import tempfile
import textwrap
from pathlib import Path

import biom
import pandas as pd

from q2_deseq2._run_data import write_run_result_artifact, DESeq2RunResult
from q2_deseq2.types import DESeq2RunDirectoryFormat


def _prepare_inputs(
    table: biom.Table,
    condition,
    min_total_count: int,
    test_level: str,
    reference_level: str,
) -> tuple[pd.DataFrame, pd.DataFrame, str, str]:
    counts = table.to_dataframe(dense=True)
    sample_ids = list(table.ids(axis="sample"))

    metadata = condition.to_series().dropna().astype(str)
    matched_samples = [
        sample_id for sample_id in sample_ids if sample_id in metadata.index
    ]
    if len(matched_samples) < 2:
        raise ValueError(
            "At least two samples must overlap between the feature table and condition metadata."
        )

    counts = counts.loc[:, matched_samples]
    metadata = metadata.loc[matched_samples]

    if metadata.nunique() < 2:
        raise ValueError("Condition metadata must contain at least two unique levels.")

    counts = counts.round().astype(int)
    if (counts.values < 0).any():
        raise ValueError("Feature table counts must be non-negative integers.")

    counts = counts.loc[counts.sum(axis=1) >= min_total_count]
    if counts.empty:
        raise ValueError(
            "No genes remain after filtering. Lower min_total_count or provide a denser feature table."
        )

    has_test = bool(test_level)
    has_reference = bool(reference_level)
    if has_test != has_reference:
        raise ValueError(
            "Provide both test_level and reference_level, or leave both unset."
        )

    levels = sorted(metadata.unique().tolist())
    if not has_test:
        if len(levels) != 2:
            raise ValueError(
                "Condition metadata has more than two levels. "
                "Set both test_level and reference_level to define a contrast."
            )
        reference_level = levels[0]
        test_level = levels[1]
    else:
        if test_level not in levels:
            raise ValueError(
                f'test_level "{test_level}" is not present in the condition metadata.'
            )
        if reference_level not in levels:
            raise ValueError(
                f'reference_level "{reference_level}" is not present in the condition metadata.'
            )
        if test_level == reference_level:
            raise ValueError("test_level and reference_level must be different values.")

    coldata = pd.DataFrame({"condition": metadata})
    return counts, coldata, test_level, reference_level


def _write_r_script(script_fp: Path) -> None:
    script = textwrap.dedent(
        """
        args <- commandArgs(trailingOnly = TRUE)

        get_arg <- function(flag) {
          idx <- match(flag, args)
          if (is.na(idx) || idx == length(args)) {
            stop(paste("Missing argument for", flag))
          }
          args[[idx + 1]]
        }

        parse_bool <- function(value) {
          lower <- tolower(value)
          lower %in% c("true", "1", "yes")
        }

        counts_path <- get_arg("--counts")
        coldata_path <- get_arg("--coldata")
        results_path <- get_arg("--results")
        norm_counts_path <- get_arg("--normalized-counts")
        summary_path <- get_arg("--summary")
        ma_plot_path <- get_arg("--ma-plot")
        volcano_plot_path <- get_arg("--volcano-plot")
        fit_type <- get_arg("--fit-type")
        alpha <- as.numeric(get_arg("--alpha"))
        cooks_cutoff <- parse_bool(get_arg("--cooks-cutoff"))
        independent_filtering <- parse_bool(get_arg("--independent-filtering"))
        test_level <- get_arg("--test-level")
        reference_level <- get_arg("--reference-level")

        suppressPackageStartupMessages(library("DESeq2"))

        counts <- read.table(
          counts_path,
          header = TRUE,
          sep = "\\t",
          row.names = 1,
          check.names = FALSE,
          quote = "",
          comment.char = ""
        )
        coldata <- read.table(
          coldata_path,
          header = TRUE,
          sep = "\\t",
          row.names = 1,
          check.names = FALSE,
          quote = "",
          comment.char = ""
        )

        counts <- round(as.matrix(counts))
        storage.mode(counts) <- "integer"
        coldata$condition <- factor(coldata$condition)

        if (!identical(colnames(counts), rownames(coldata))) {
          stop("Sample order mismatch between count matrix and metadata.")
        }

        dds <- DESeqDataSetFromMatrix(
          countData = counts,
          colData = coldata,
          design = ~ condition
        )
        dds <- dds[rowSums(counts(dds)) > 0, ]
        dds <- DESeq(dds, fitType = fit_type)

        if (test_level != "" && reference_level != "") {
          res <- results(
            dds,
            contrast = c("condition", test_level, reference_level),
            alpha = alpha,
            cooksCutoff = cooks_cutoff,
            independentFiltering = independent_filtering
          )
        } else {
          res <- results(
            dds,
            alpha = alpha,
            cooksCutoff = cooks_cutoff,
            independentFiltering = independent_filtering
          )
        }

        res_df <- as.data.frame(res)
        res_df$feature_id <- rownames(res_df)
        res_df <- res_df[, c("feature_id", setdiff(colnames(res_df), "feature_id"))]

        if ("padj" %in% colnames(res_df) && "pvalue" %in% colnames(res_df)) {
          res_df <- res_df[order(is.na(res_df$padj), res_df$padj, is.na(res_df$pvalue), res_df$pvalue), ]
        }

        write.table(
          res_df,
          file = results_path,
          sep = "\\t",
          quote = FALSE,
          row.names = FALSE
        )

        norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
        norm_counts$feature_id <- rownames(norm_counts)
        norm_counts <- norm_counts[, c("feature_id", setdiff(colnames(norm_counts), "feature_id"))]
        write.table(
          norm_counts,
          file = norm_counts_path,
          sep = "\\t",
          quote = FALSE,
          row.names = FALSE
        )

        writeLines(capture.output(summary(res)), con = summary_path)

        plot_label <- if (test_level != "" && reference_level != "") {
          sprintf("%s vs %s", test_level, reference_level)
        } else {
          "default contrast"
        }
        png(filename = ma_plot_path, width = 1200, height = 900)
        plotMA(res, alpha = alpha, main = paste("DESeq2 MA plot:", plot_label))
        dev.off()

        y_values <- -log10(res_df$padj)
        finite_idx <- is.finite(res_df$log2FoldChange) & is.finite(y_values)

        png(filename = volcano_plot_path, width = 1200, height = 900)
        if (any(finite_idx)) {
          plot(
            res_df$log2FoldChange[finite_idx],
            y_values[finite_idx],
            pch = 20,
            col = rgb(0.2, 0.4, 0.7, 0.65),
            xlab = "log2 fold change",
            ylab = "-log10 adjusted p-value",
            main = paste("DESeq2 Volcano Plot:", plot_label)
          )
          abline(v = 0, col = "#1F77B4", lty = 2)
          if (alpha > 0 && alpha < 1) {
            abline(h = -log10(alpha), col = "#D62728", lty = 3)
          }
        } else {
          plot.new()
          title(main = paste("DESeq2 Volcano Plot:", plot_label))
          text(0.5, 0.5, "No finite points available for volcano plot.")
        }
        dev.off()
        """
    ).strip()
    script_fp.write_text(script + "\n", encoding="utf-8")


def run_deseq2(
    table: biom.Table,
    condition,
    test_level: str = "",
    reference_level: str = "",
    min_total_count: int = 10,
    fit_type: str = "parametric",
    alpha: float = 0.05,
    cooks_cutoff: bool = True,
    independent_filtering: bool = True,
) -> DESeq2RunResult:
    counts_df, coldata_df, test_level, reference_level = _prepare_inputs(
        table, condition, min_total_count, test_level, reference_level
    )

    with tempfile.TemporaryDirectory(prefix="q2-deseq2-") as temp_dir:
        temp_path = Path(temp_dir)
        counts_fp = temp_path / "counts.tsv"
        coldata_fp = temp_path / "coldata.tsv"
        results_fp = temp_path / "deseq2_results.tsv"
        normalized_counts_fp = temp_path / "normalized_counts.tsv"
        summary_fp = temp_path / "deseq2_summary.txt"
        ma_plot_fp = temp_path / "ma_plot.png"
        volcano_plot_fp = temp_path / "volcano_plot.png"
        script_fp = temp_path / "run_deseq2.R"

        counts_df.to_csv(counts_fp, sep="\t", index_label="feature_id")
        coldata_df.to_csv(coldata_fp, sep="\t", index_label="sample_id")
        _write_r_script(script_fp)

        cmd = [
            "Rscript",
            str(script_fp),
            "--counts",
            str(counts_fp),
            "--coldata",
            str(coldata_fp),
            "--results",
            str(results_fp),
            "--normalized-counts",
            str(normalized_counts_fp),
            "--summary",
            str(summary_fp),
            "--ma-plot",
            str(ma_plot_fp),
            "--volcano-plot",
            str(volcano_plot_fp),
            "--fit-type",
            fit_type,
            "--alpha",
            str(alpha),
            "--cooks-cutoff",
            str(cooks_cutoff).lower(),
            "--independent-filtering",
            str(independent_filtering).lower(),
            "--test-level",
            test_level,
            "--reference-level",
            reference_level,
        ]

        try:
            run(cmd, check=True, capture_output=True, text=True)
        except CalledProcessError as exc:
            detail = exc.stderr.strip() or exc.stdout.strip()
            raise RuntimeError(
                f"DESeq2 command failed with exit code {exc.returncode}: {detail}"
            ) from exc

        expected_outputs = {
            "deseq2_results.tsv": results_fp,
            "normalized_counts.tsv": normalized_counts_fp,
            "ma_plot.png": ma_plot_fp,
            "volcano_plot.png": volcano_plot_fp,
        }
        for expected_name, path in expected_outputs.items():
            if not path.exists():
                raise RuntimeError(
                    f'DESeq2 completed but expected output file "{expected_name}" was not created.'
                )

        results_df = pd.read_csv(results_fp, sep="\t")
        normalized_counts_df = pd.read_csv(normalized_counts_fp, sep="\t")
        ma_plot_png = ma_plot_fp.read_bytes()
        volcano_plot_png = volcano_plot_fp.read_bytes()

        return DESeq2RunResult(
            results=results_df,
            normalized_counts=normalized_counts_df,
            ma_plot_png=ma_plot_png,
            volcano_plot_png=volcano_plot_png,
            test_level=test_level,
            reference_level=reference_level,
        )


def _estimate_differential_expression(
    table: biom.Table,
    condition: str,
    test_level: str = "",
    reference_level: str = "",
    min_total_count: int = 10,
    fit_type: str = "parametric",
    alpha: float = 0.05,
    cooks_cutoff: bool = True,
    independent_filtering: bool = True,
) -> (pd.DataFrame, DESeq2RunDirectoryFormat):
    run_result = run_deseq2(
        table=table,
        condition=condition,
        test_level=test_level,
        reference_level=reference_level,
        min_total_count=min_total_count,
        fit_type=fit_type,
        alpha=alpha,
        cooks_cutoff=cooks_cutoff,
        independent_filtering=independent_filtering,
    )
    run_data = write_run_result_artifact(run_result=run_result, alpha=alpha)
    return run_result.results, run_data
