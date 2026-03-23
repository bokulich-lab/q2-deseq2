# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
from pathlib import Path
from subprocess import CalledProcessError, run
import tempfile
import textwrap

import biom
import pandas as pd
from pandas.api.types import is_numeric_dtype

from q2_deseq2._run_data import DESeq2RunResult, write_run_result_artifact
from q2_deseq2.types import DESeq2RunDirectoryFormat

_FORMULA_ALLOWED_CHARS = re.compile(r"^[A-Za-z0-9_:+*() \t]+$")
_FORMULA_TOKEN_RE = re.compile(r"\b[A-Za-z_][A-Za-z0-9_]*\b")
_COEF_SPEC_RE = re.compile(r"^coef::(.+)$")
_CONTRAST_SPEC_RE = re.compile(r"^contrast::([^:]+)::([^:]+)::([^:]+)$")
_SIMPLE_SPEC_RE = re.compile(
    r"^simple::([^:]+)::([^:]+)::([^:]+)\|within::([^:]+)::([^:]+)$"
)


def _prepare_count_table(table: biom.Table) -> pd.DataFrame:
    counts = table.to_dataframe(dense=True).round().astype(int)
    if (counts.values < 0).any():
        raise ValueError("Feature table counts must be non-negative integers.")
    return counts


def _filter_counts(counts: pd.DataFrame, min_total_count: int) -> pd.DataFrame:
    filtered = counts.loc[counts.sum(axis=1) >= min_total_count]
    if filtered.empty:
        raise ValueError(
            "No genes remain after filtering. Lower min_total_count or provide a denser feature table."
        )
    return filtered


def _collect_matching_samples(
    sample_ids: list[str], metadata_index: pd.Index, context: str
) -> list[str]:
    matched_samples = [sample_id for sample_id in sample_ids if sample_id in metadata_index]
    if len(matched_samples) < 2:
        raise ValueError(
            f"At least two samples must overlap between the feature table and {context}."
        )
    return matched_samples


def _normalize_formula(formula: str, parameter_name: str) -> str:
    normalized = str(formula).strip()
    if not normalized:
        raise ValueError(f"{parameter_name} is required.")
    if not _FORMULA_ALLOWED_CHARS.fullmatch(normalized):
        raise ValueError(
            f"{parameter_name} may only contain metadata column names, spaces, '+', ':', '*', and parentheses."
        )
    return normalized


def _extract_formula_columns(formula: str, parameter_name: str) -> list[str]:
    normalized = _normalize_formula(formula, parameter_name)
    columns = list(dict.fromkeys(_FORMULA_TOKEN_RE.findall(normalized)))
    if not columns:
        raise ValueError(
            f"{parameter_name} must reference at least one metadata column."
        )
    return columns


def _parse_reference_level_spec(reference_level_spec: str) -> tuple[str, str]:
    column, separator, level = str(reference_level_spec).strip().partition("::")
    if separator != "::" or not column or not level:
        raise ValueError(
            "reference_levels entries must use the form 'column::level'."
        )
    return column, level


def _normalize_reference_levels(reference_levels: list[str] | None) -> list[str]:
    normalized = []
    seen_columns = set()
    for raw_spec in reference_levels or []:
        spec = str(raw_spec).strip()
        if not spec:
            continue
        column, level = _parse_reference_level_spec(spec)
        if column in seen_columns:
            raise ValueError(
                f'Duplicate reference_levels entry provided for metadata column "{column}".'
            )
        normalized.append(f"{column}::{level}")
        seen_columns.add(column)
    return normalized


def _validate_effect_specs(effect_specs: list[str] | None, test: str) -> list[str]:
    normalized = []
    for raw_spec in effect_specs or []:
        spec = str(raw_spec).strip()
        if not spec:
            continue
        if (
            _COEF_SPEC_RE.match(spec) is None
            and _CONTRAST_SPEC_RE.match(spec) is None
            and _SIMPLE_SPEC_RE.match(spec) is None
        ):
            raise ValueError(
                "effect_specs entries must use one of: "
                "'coef::<resultsName>', "
                "'contrast::<factor>::<numerator>::<denominator>', or "
                "'simple::<factor>::<numerator>::<denominator>|within::<factor>::<level>'."
            )
        normalized.append(spec)

    if test == "lrt" and normalized:
        raise ValueError("effect_specs are not supported when test='lrt'.")

    return normalized


def _coerce_metadata_column(column: pd.Series) -> pd.Series:
    if is_numeric_dtype(column):
        return pd.to_numeric(column, errors="raise")
    return column.astype(str)


def _prepare_inputs(
    table: biom.Table,
    condition,
    min_total_count: int,
    reference_level: str,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], str]:
    counts = _prepare_count_table(table)
    sample_ids = list(table.ids(axis="sample"))

    metadata = condition.to_series().dropna().astype(str)
    matched_samples = _collect_matching_samples(
        sample_ids, metadata.index, "condition metadata"
    )

    counts = counts.loc[:, matched_samples]
    metadata = metadata.loc[matched_samples]

    if metadata.nunique() < 2:
        raise ValueError("Condition metadata must contain at least two unique levels.")

    counts = _filter_counts(counts, min_total_count)

    levels = sorted(metadata.unique().tolist())
    if not reference_level:
        if len(levels) != 2:
            raise ValueError(
                "Condition metadata has more than two levels. "
                "Set reference_level to define the baseline for all pairwise contrasts."
            )
        reference_level = levels[0]
    elif reference_level not in levels:
        raise ValueError(
            f'reference_level "{reference_level}" is not present in the condition metadata.'
        )

    comparison_levels = [level for level in levels if level != reference_level]
    if not comparison_levels:
        raise ValueError(
            "Condition metadata must contain at least one non-reference level."
        )

    coldata = pd.DataFrame({"condition": metadata})
    return counts, coldata, comparison_levels, reference_level


def _prepare_model_inputs(
    table: biom.Table,
    metadata,
    fixed_effects_formula: str,
    min_total_count: int,
    reference_levels: list[str] | None = None,
    reduced_formula: str = "",
) -> tuple[pd.DataFrame, pd.DataFrame, str, list[str], str]:
    normalized_formula = _normalize_formula(
        fixed_effects_formula, parameter_name="fixed_effects_formula"
    )
    normalized_reduced_formula = (
        _normalize_formula(reduced_formula, parameter_name="reduced_formula")
        if str(reduced_formula).strip()
        else ""
    )

    metadata_df = metadata.to_dataframe()
    sample_ids = list(table.ids(axis="sample"))
    counts = _prepare_count_table(table)
    referenced_columns = _extract_formula_columns(
        normalized_formula, parameter_name="fixed_effects_formula"
    )
    if normalized_reduced_formula:
        referenced_columns.extend(
            _extract_formula_columns(
                normalized_reduced_formula, parameter_name="reduced_formula"
            )
        )
    referenced_columns = list(dict.fromkeys(referenced_columns))

    missing_columns = [
        column for column in referenced_columns if column not in metadata_df.columns
    ]
    if missing_columns:
        missing_display = ", ".join(sorted(missing_columns))
        raise ValueError(
            "The metadata is missing columns required by the model formula: "
            f"{missing_display}."
        )

    normalized_reference_levels = _normalize_reference_levels(reference_levels)
    for reference_level_spec in normalized_reference_levels:
        column, _ = _parse_reference_level_spec(reference_level_spec)
        if column not in referenced_columns:
            raise ValueError(
                f'reference_levels entry "{reference_level_spec}" refers to "{column}", '
                "but that column is not used in the fitted model."
            )

    matched_samples = _collect_matching_samples(
        sample_ids, metadata_df.index, "metadata"
    )
    coldata = metadata_df.loc[matched_samples, referenced_columns].copy()
    counts = counts.loc[:, matched_samples]

    coldata = coldata.loc[~coldata.isna().any(axis=1)].copy()
    if coldata.shape[0] < 2:
        raise ValueError(
            "At least two samples with complete metadata are required after dropping missing values."
        )

    counts = counts.loc[:, coldata.index]
    counts = _filter_counts(counts, min_total_count)
    for column in coldata.columns:
        coldata[column] = _coerce_metadata_column(coldata[column])

    for reference_level_spec in normalized_reference_levels:
        column, level = _parse_reference_level_spec(reference_level_spec)
        if is_numeric_dtype(coldata[column]):
            raise ValueError(
                f'reference_levels entry "{reference_level_spec}" refers to numeric metadata column "{column}".'
            )
        observed_levels = set(coldata[column].astype(str))
        if level not in observed_levels:
            raise ValueError(
                f'reference_levels entry "{reference_level_spec}" refers to level "{level}", '
                f'which is not present in metadata column "{column}".'
            )

    return (
        counts,
        coldata,
        normalized_formula,
        normalized_reference_levels,
        normalized_reduced_formula,
    )


def _write_lines(path: Path, lines: list[str]) -> None:
    payload = "\n".join(lines)
    if payload:
        payload += "\n"
    path.write_text(payload, encoding="utf-8")


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

        read_list_file <- function(path) {
          if (!file.exists(path)) {
            return(character())
          }
          values <- readLines(path, warn = FALSE)
          values[nzchar(trimws(values))]
        }

        parse_reference_level <- function(spec) {
          parts <- strsplit(spec, "::", fixed = TRUE)[[1]]
          if (length(parts) != 2 || any(parts == "")) {
            stop(sprintf("Invalid reference_levels entry: %s", spec))
          }
          parts
        }

        ordered_levels <- function(values, reference_level) {
          all_levels <- sort(unique(as.character(values)))
          if (!(reference_level %in% all_levels)) {
            stop(sprintf('Reference level "%s" is not present.', reference_level))
          }
          c(reference_level, setdiff(all_levels, reference_level))
        }

        coerce_coldata <- function(coldata) {
          for (column_name in names(coldata)) {
            if (is.character(coldata[[column_name]])) {
              coldata[[column_name]] <- factor(coldata[[column_name]])
            }
          }
          coldata
        }

        apply_reference_levels <- function(coldata, specs) {
          if (length(specs) == 0) {
            return(coerce_coldata(coldata))
          }

          coldata <- coerce_coldata(coldata)
          for (spec in specs) {
            parsed <- parse_reference_level(spec)
            column_name <- parsed[[1]]
            reference_level <- parsed[[2]]

            if (!(column_name %in% names(coldata))) {
              stop(sprintf('Reference level column "%s" is not present in the metadata.', column_name))
            }
            if (is.numeric(coldata[[column_name]])) {
              stop(sprintf('Reference level column "%s" is numeric and cannot be releveled.', column_name))
            }

            values <- as.character(coldata[[column_name]])
            coldata[[column_name]] <- factor(
              values,
              levels = ordered_levels(values, reference_level)
            )
          }

          coldata
        }

        relevel_for_simple_effect <- function(coldata, factor_name, denominator, within_factor, within_level) {
          if (factor_name == within_factor) {
            stop("simple effect target factor and within factor must be different.")
          }

          for (column_name in c(factor_name, within_factor)) {
            if (!(column_name %in% names(coldata))) {
              stop(sprintf('Metadata column "%s" required for simple effect is not present.', column_name))
            }
            if (is.numeric(coldata[[column_name]])) {
              stop(sprintf('Metadata column "%s" used in a simple effect must be categorical.', column_name))
            }
          }

          releveled <- coldata
          factor_values <- as.character(releveled[[factor_name]])
          within_values <- as.character(releveled[[within_factor]])

          releveled[[factor_name]] <- factor(
            factor_values,
            levels = ordered_levels(factor_values, denominator)
          )
          releveled[[within_factor]] <- factor(
            within_values,
            levels = ordered_levels(within_values, within_level)
          )
          releveled
        }

        make_result_frame <- function(res, effect_id, effect_label, effect_kind, effect_expression, comparison = "", test_level = "", reference_level = "") {
          res_df <- as.data.frame(res)
          res_df$effect_id <- effect_id
          res_df$effect_label <- effect_label
          res_df$effect_kind <- effect_kind
          res_df$effect_expression <- effect_expression
          res_df$comparison <- comparison
          res_df$test_level <- test_level
          res_df$reference_level <- reference_level
          res_df$feature_id <- rownames(res_df)
          preferred_columns <- c(
            "effect_id",
            "effect_label",
            "effect_kind",
            "effect_expression",
            "comparison",
            "test_level",
            "reference_level",
            "feature_id"
          )
          res_df <- res_df[
            ,
            c(preferred_columns, setdiff(colnames(res_df), preferred_columns))
          ]

          if ("padj" %in% colnames(res_df) && "pvalue" %in% colnames(res_df)) {
            res_df <- res_df[
              order(is.na(res_df$padj), res_df$padj, is.na(res_df$pvalue), res_df$pvalue),
            ]
          }

          res_df
        }

        counts_path <- get_arg("--counts")
        coldata_path <- get_arg("--coldata")
        results_path <- get_arg("--results")
        norm_counts_path <- get_arg("--normalized-counts")
        summary_path <- get_arg("--summary")
        ma_plot_path <- get_arg("--ma-plot")
        volcano_plot_path <- get_arg("--volcano-plot")
        results_names_path <- get_arg("--results-names")
        reference_levels_path <- get_arg("--reference-levels")
        effect_specs_path <- get_arg("--effect-specs")
        fit_type <- get_arg("--fit-type")
        alpha <- as.numeric(get_arg("--alpha"))
        cooks_cutoff <- parse_bool(get_arg("--cooks-cutoff"))
        independent_filtering <- parse_bool(get_arg("--independent-filtering"))
        fixed_effects_formula <- get_arg("--fixed-effects-formula")
        test_kind <- tolower(get_arg("--test"))
        reduced_formula <- get_arg("--reduced-formula")

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
        coldata <- apply_reference_levels(coldata, read_list_file(reference_levels_path))

        if (!identical(colnames(counts), rownames(coldata))) {
          stop("Sample order mismatch between count matrix and metadata.")
        }

        design_formula <- as.formula(paste("~", fixed_effects_formula))
        dds <- DESeqDataSetFromMatrix(
          countData = counts,
          colData = coldata,
          design = design_formula
        )
        dds <- dds[rowSums(counts(dds)) > 0, ]

        if (test_kind == "lrt") {
          if (!nzchar(trimws(reduced_formula))) {
            stop("reduced_formula is required when test='lrt'.")
          }
          if (length(read_list_file(effect_specs_path)) > 0) {
            stop("effect_specs are not supported when test='lrt'.")
          }
          dds <- DESeq(
            dds,
            test = "LRT",
            reduced = as.formula(paste("~", reduced_formula)),
            fitType = fit_type
          )
        } else {
          dds <- DESeq(dds, fitType = fit_type)
        }

        results_names <- resultsNames(dds)
        writeLines(results_names, con = results_names_path)

        result_frames <- list()
        summary_lines <- character()
        default_res <- NULL
        default_res_df <- NULL
        default_effect_label <- NULL

        if (test_kind == "lrt") {
          effect_id <- sprintf("lrt::%s::vs::%s", fixed_effects_formula, reduced_formula)
          effect_label <- sprintf("LRT: %s vs %s", fixed_effects_formula, reduced_formula)
          res <- results(
            dds,
            alpha = alpha,
            cooksCutoff = cooks_cutoff,
            independentFiltering = independent_filtering
          )
          res_df <- make_result_frame(
            res,
            effect_id = effect_id,
            effect_label = effect_label,
            effect_kind = "lrt",
            effect_expression = sprintf("full=~ %s; reduced=~ %s", fixed_effects_formula, reduced_formula)
          )
          result_frames[[effect_id]] <- res_df
          summary_lines <- c(
            summary_lines,
            effect_label,
            capture.output(summary(res)),
            ""
          )
          default_res <- res
          default_res_df <- res_df
          default_effect_label <- effect_label
        } else {
          effect_specs <- read_list_file(effect_specs_path)
          if (length(effect_specs) == 0) {
            coefficient_names <- setdiff(results_names, "Intercept")
            effect_specs <- sprintf("coef::%s", coefficient_names)
          }
          if (length(effect_specs) == 0) {
            stop("No non-intercept coefficients were available to extract from the fitted model.")
          }

          simple_dds_cache <- list()

          for (spec in effect_specs) {
            if (startsWith(spec, "coef::")) {
              coefficient_name <- substring(spec, nchar("coef::") + 1)
              if (!(coefficient_name %in% results_names)) {
                stop(sprintf('Requested coefficient "%s" is not available in resultsNames(dds).', coefficient_name))
              }
              res <- results(
                dds,
                name = coefficient_name,
                alpha = alpha,
                cooksCutoff = cooks_cutoff,
                independentFiltering = independent_filtering
              )
              effect_label <- coefficient_name
              res_df <- make_result_frame(
                res,
                effect_id = spec,
                effect_label = effect_label,
                effect_kind = "coefficient",
                effect_expression = sprintf('name="%s"', coefficient_name)
              )
            } else if (startsWith(spec, "contrast::")) {
              parsed <- strsplit(spec, "::", fixed = TRUE)[[1]]
              factor_name <- parsed[[2]]
              numerator <- parsed[[3]]
              denominator <- parsed[[4]]
              if (!(factor_name %in% names(coldata))) {
                stop(sprintf('Contrast factor "%s" is not present in the metadata.', factor_name))
              }
              if (is.numeric(coldata[[factor_name]])) {
                stop(sprintf('Contrast factor "%s" must be categorical.', factor_name))
              }
              res <- results(
                dds,
                contrast = c(factor_name, numerator, denominator),
                alpha = alpha,
                cooksCutoff = cooks_cutoff,
                independentFiltering = independent_filtering
              )
              comparison_label <- sprintf("%s vs. %s", numerator, denominator)
              effect_label <- sprintf("%s: %s vs %s", factor_name, numerator, denominator)
              res_df <- make_result_frame(
                res,
                effect_id = spec,
                effect_label = effect_label,
                effect_kind = "contrast",
                effect_expression = sprintf('contrast=c("%s","%s","%s")', factor_name, numerator, denominator),
                comparison = comparison_label,
                test_level = numerator,
                reference_level = denominator
              )
            } else if (startsWith(spec, "simple::")) {
              parsed <- regmatches(
                spec,
                regexec("^simple::([^:]+)::([^:]+)::([^:]+)\\\\|within::([^:]+)::([^:]+)$", spec)
              )[[1]]
              if (length(parsed) != 6) {
                stop(sprintf("Invalid simple effect spec: %s", spec))
              }
              factor_name <- parsed[[2]]
              numerator <- parsed[[3]]
              denominator <- parsed[[4]]
              within_factor <- parsed[[5]]
              within_level <- parsed[[6]]
              cache_key <- sprintf("%s::%s::%s::%s", factor_name, denominator, within_factor, within_level)

              if (is.null(simple_dds_cache[[cache_key]])) {
                simple_coldata <- relevel_for_simple_effect(
                  coldata,
                  factor_name = factor_name,
                  denominator = denominator,
                  within_factor = within_factor,
                  within_level = within_level
                )
                simple_dds <- DESeqDataSetFromMatrix(
                  countData = counts,
                  colData = simple_coldata,
                  design = design_formula
                )
                simple_dds <- simple_dds[rowSums(counts(simple_dds)) > 0, ]
                simple_dds <- DESeq(simple_dds, fitType = fit_type)
                simple_dds_cache[[cache_key]] <- simple_dds
              }

              simple_dds <- simple_dds_cache[[cache_key]]
              res <- results(
                simple_dds,
                contrast = c(factor_name, numerator, denominator),
                alpha = alpha,
                cooksCutoff = cooks_cutoff,
                independentFiltering = independent_filtering
              )
              comparison_label <- sprintf("%s vs. %s (%s=%s)", numerator, denominator, within_factor, within_level)
              effect_label <- sprintf(
                "%s: %s vs %s within %s=%s",
                factor_name,
                numerator,
                denominator,
                within_factor,
                within_level
              )
              res_df <- make_result_frame(
                res,
                effect_id = spec,
                effect_label = effect_label,
                effect_kind = "simple",
                effect_expression = sprintf(
                  'contrast=c("%s","%s","%s"); within="%s=%s"',
                  factor_name,
                  numerator,
                  denominator,
                  within_factor,
                  within_level
                ),
                comparison = comparison_label,
                test_level = numerator,
                reference_level = denominator
              )
            } else {
              stop(sprintf("Unsupported effect spec: %s", spec))
            }

            result_frames[[spec]] <- res_df
            summary_lines <- c(
              summary_lines,
              res_df$effect_label[[1]],
              capture.output(summary(res)),
              ""
            )

            if (is.null(default_res)) {
              default_res <- res
              default_res_df <- res_df
              default_effect_label <- res_df$effect_label[[1]]
            }
          }
        }

        res_df <- do.call(rbind, result_frames)

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

        if (length(summary_lines) == 0) {
          summary_lines <- "No effects were generated."
        }
        writeLines(summary_lines, con = summary_path)

        plot_label <- default_effect_label
        png(filename = ma_plot_path, width = 1200, height = 900)
        plotMA(default_res, alpha = alpha, main = paste("DESeq2 MA plot:", plot_label))
        dev.off()

        y_values <- -log10(default_res_df$padj)
        finite_idx <- is.finite(default_res_df$log2FoldChange) & is.finite(y_values)

        png(filename = volcano_plot_path, width = 1200, height = 900)
        if (any(finite_idx)) {
          plot(
            default_res_df$log2FoldChange[finite_idx],
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


def _first_non_empty_string(value) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return text


def _first_value_from_column(frame: pd.DataFrame, column_name: str) -> str:
    if column_name not in frame.columns:
        return ""
    for value in frame[column_name]:
        text = _first_non_empty_string(value)
        if text:
            return text
    return ""


def _unique_non_empty_values(values) -> tuple[str, ...]:
    normalized = []
    seen = set()
    for value in values:
        text = _first_non_empty_string(value)
        if text and text not in seen:
            normalized.append(text)
            seen.add(text)
    return tuple(normalized)


def _run_deseq2_with_frames(
    counts_df: pd.DataFrame,
    coldata_df: pd.DataFrame,
    fixed_effects_formula: str,
    reference_levels: list[str] | None = None,
    effect_specs: list[str] | None = None,
    test: str = "wald",
    reduced_formula: str = "",
    fit_type: str = "parametric",
    alpha: float = 0.05,
    cooks_cutoff: bool = True,
    independent_filtering: bool = True,
    legacy_test_level: str = "",
    legacy_reference_level: str = "",
) -> DESeq2RunResult:
    normalized_reference_levels = _normalize_reference_levels(reference_levels)
    normalized_test = str(test).strip().lower() or "wald"
    if normalized_test not in {"wald", "lrt"}:
        raise ValueError("test must be either 'wald' or 'lrt'.")
    normalized_effect_specs = _validate_effect_specs(effect_specs, normalized_test)
    normalized_reduced_formula = (
        _normalize_formula(reduced_formula, parameter_name="reduced_formula")
        if str(reduced_formula).strip()
        else ""
    )
    if normalized_test == "lrt" and not normalized_reduced_formula:
        raise ValueError("reduced_formula is required when test='lrt'.")

    counts_df = counts_df + 1

    with tempfile.TemporaryDirectory(prefix="q2-deseq2-") as temp_dir:
        temp_path = Path(temp_dir)
        counts_fp = temp_path / "counts.tsv"
        coldata_fp = temp_path / "coldata.tsv"
        results_fp = temp_path / "deseq2_results.tsv"
        normalized_counts_fp = temp_path / "normalized_counts.tsv"
        summary_fp = temp_path / "deseq2_summary.txt"
        ma_plot_fp = temp_path / "ma_plot.png"
        volcano_plot_fp = temp_path / "volcano_plot.png"
        results_names_fp = temp_path / "results_names.txt"
        reference_levels_fp = temp_path / "reference_levels.txt"
        effect_specs_fp = temp_path / "effect_specs.txt"
        script_fp = temp_path / "run_deseq2.R"

        counts_df.to_csv(counts_fp, sep="\t", index_label="feature_id")
        coldata_df.to_csv(coldata_fp, sep="\t", index_label="sample_id")
        _write_lines(reference_levels_fp, normalized_reference_levels)
        _write_lines(effect_specs_fp, normalized_effect_specs)
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
            "--results-names",
            str(results_names_fp),
            "--reference-levels",
            str(reference_levels_fp),
            "--effect-specs",
            str(effect_specs_fp),
            "--fit-type",
            fit_type,
            "--alpha",
            str(alpha),
            "--cooks-cutoff",
            str(cooks_cutoff).lower(),
            "--independent-filtering",
            str(independent_filtering).lower(),
            "--fixed-effects-formula",
            fixed_effects_formula,
            "--test",
            normalized_test,
            "--reduced-formula",
            normalized_reduced_formula,
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
            "results_names.txt": results_names_fp,
        }
        for expected_name, path in expected_outputs.items():
            if not path.exists():
                raise RuntimeError(
                    f'DESeq2 completed but expected output file "{expected_name}" was not created.'
                )

        results_df = pd.read_csv(results_fp, sep="\t")
        normalized_counts_df = pd.read_csv(normalized_counts_fp, sep="\t")
        available_results_names = _unique_non_empty_values(
            results_names_fp.read_text(encoding="utf-8").splitlines()
        )
        selected_effect_specs = _unique_non_empty_values(
            results_df["effect_id"] if "effect_id" in results_df.columns else ()
        )
        default_effect_id = _first_value_from_column(results_df, "effect_id")

        return DESeq2RunResult(
            results=results_df,
            normalized_counts=normalized_counts_df,
            ma_plot_png=ma_plot_fp.read_bytes(),
            volcano_plot_png=volcano_plot_fp.read_bytes(),
            test_level=legacy_test_level or _first_value_from_column(results_df, "test_level"),
            reference_level=legacy_reference_level
            or _first_value_from_column(results_df, "reference_level"),
            default_effect_id=default_effect_id,
            fixed_effects_formula=fixed_effects_formula,
            reference_levels=tuple(normalized_reference_levels),
            test=normalized_test,
            reduced_formula=normalized_reduced_formula,
            available_results_names=available_results_names,
            selected_effect_specs=selected_effect_specs,
        )


def run_deseq2_model(
    table: biom.Table,
    metadata,
    fixed_effects_formula: str,
    reference_levels: list[str] | None = None,
    effect_specs: list[str] | None = None,
    test: str = "wald",
    reduced_formula: str = "",
    min_total_count: int = 10,
    fit_type: str = "parametric",
    alpha: float = 0.05,
    cooks_cutoff: bool = True,
    independent_filtering: bool = True,
) -> DESeq2RunResult:
    (
        counts_df,
        coldata_df,
        normalized_formula,
        normalized_reference_levels,
        normalized_reduced_formula,
    ) = _prepare_model_inputs(
        table=table,
        metadata=metadata,
        fixed_effects_formula=fixed_effects_formula,
        min_total_count=min_total_count,
        reference_levels=reference_levels,
        reduced_formula=reduced_formula,
    )

    return _run_deseq2_with_frames(
        counts_df=counts_df,
        coldata_df=coldata_df,
        fixed_effects_formula=normalized_formula,
        reference_levels=normalized_reference_levels,
        effect_specs=effect_specs,
        test=test,
        reduced_formula=normalized_reduced_formula,
        fit_type=fit_type,
        alpha=alpha,
        cooks_cutoff=cooks_cutoff,
        independent_filtering=independent_filtering,
    )


def run_deseq2(
    table: biom.Table,
    condition,
    reference_level: str = "",
    min_total_count: int = 10,
    fit_type: str = "parametric",
    alpha: float = 0.05,
    cooks_cutoff: bool = True,
    independent_filtering: bool = True,
) -> DESeq2RunResult:
    counts_df, coldata_df, comparison_levels, reference_level = _prepare_inputs(
        table, condition, min_total_count, reference_level
    )
    effect_specs = [
        f"contrast::condition::{comparison_level}::{reference_level}"
        for comparison_level in comparison_levels
    ]

    return _run_deseq2_with_frames(
        counts_df=counts_df,
        coldata_df=coldata_df,
        fixed_effects_formula="condition",
        reference_levels=[f"condition::{reference_level}"],
        effect_specs=effect_specs,
        test="wald",
        fit_type=fit_type,
        alpha=alpha,
        cooks_cutoff=cooks_cutoff,
        independent_filtering=independent_filtering,
        legacy_test_level=comparison_levels[0],
        legacy_reference_level=reference_level,
    )


def _estimate_model(
    table: biom.Table,
    metadata: str,
    fixed_effects_formula: str,
    reference_levels: list[str] | None = None,
    effect_specs: list[str] | None = None,
    test: str = "wald",
    reduced_formula: str = "",
    min_total_count: int = 10,
    fit_type: str = "parametric",
    alpha: float = 0.05,
    cooks_cutoff: bool = True,
    independent_filtering: bool = True,
) -> (pd.DataFrame, DESeq2RunDirectoryFormat):
    run_result = run_deseq2_model(
        table=table,
        metadata=metadata,
        fixed_effects_formula=fixed_effects_formula,
        reference_levels=reference_levels,
        effect_specs=effect_specs,
        test=test,
        reduced_formula=reduced_formula,
        min_total_count=min_total_count,
        fit_type=fit_type,
        alpha=alpha,
        cooks_cutoff=cooks_cutoff,
        independent_filtering=independent_filtering,
    )
    run_data = write_run_result_artifact(run_result=run_result, alpha=alpha)
    return run_result.results, run_data


def _estimate_differential_expression(
    table: biom.Table,
    condition: str,
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
        reference_level=reference_level,
        min_total_count=min_total_count,
        fit_type=fit_type,
        alpha=alpha,
        cooks_cutoff=cooks_cutoff,
        independent_filtering=independent_filtering,
    )
    run_data = write_run_result_artifact(run_result=run_result, alpha=alpha)
    return run_result.results, run_data
