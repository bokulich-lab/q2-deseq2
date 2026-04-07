# ============================================================================
# run_deseq2.R — DESeq2 runner invoked by q2-deseq2 via Rscript
#
# This script is called as a subprocess by _run_deseq2_with_frames() in
# _methodlib/runner.py. It reads count and metadata TSVs written by the Python side,
# fits a DESeq2 model, extracts results for each requested effect, and writes
# output files that the Python side reads back.
# ============================================================================

# --- CLI argument helpers ---------------------------------------------------

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

# Read a file where each non-empty line is one entry (used for reference
# levels and effect specs).
read_list_file <- function(path) {
  if (!file.exists(path)) {
    return(character())
  }
  values <- readLines(path, warn = FALSE)
  values[nzchar(trimws(values))]
}

# --- Reference-level handling -----------------------------------------------

# Parse a "column::level" spec into a two-element character vector.
parse_reference_level <- function(spec) {
  parts <- strsplit(spec, "::", fixed = TRUE)[[1]]
  if (length(parts) != 2 || any(parts == "")) {
    stop(sprintf("Invalid reference_levels entry: %s", spec))
  }
  parts
}

# Return factor levels with `reference_level` first, remainder sorted.
ordered_levels <- function(values, reference_level) {
  all_levels <- sort(unique(as.character(values)))
  if (!(reference_level %in% all_levels)) {
    stop(sprintf('Reference level "%s" is not present.', reference_level))
  }
  c(reference_level, setdiff(all_levels, reference_level))
}

# Convert all character columns to factors (DESeq2 requires factors).
coerce_coldata <- function(coldata) {
  for (column_name in names(coldata)) {
    if (is.character(coldata[[column_name]])) {
      coldata[[column_name]] <- factor(coldata[[column_name]])
    }
  }
  coldata
}

# Apply user-specified reference levels by re-ordering factor levels so
# that the reference comes first (DESeq2 treats the first level as baseline).
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

# --- Simple-effect helpers --------------------------------------------------

# For "simple" effects (e.g. genotype effect within a specific treatment
# level), relevel both the target factor and the within-factor so that
# DESeq2's contrast extraction works correctly.
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

# --- Result formatting ------------------------------------------------------

# Convert a DESeq2 results object into a data.frame with metadata columns
# prepended (effect_id, effect_label, etc.) and rows sorted by padj/pvalue.
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

  # Reorder columns: metadata columns first, then DESeq2 statistics.
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

  # Sort: significant results first, then by raw p-value; NAs last.
  if ("padj" %in% colnames(res_df) && "pvalue" %in% colnames(res_df)) {
    res_df <- res_df[
      order(is.na(res_df$padj), res_df$padj, is.na(res_df$pvalue), res_df$pvalue),
    ]
  }

  res_df
}

# Run DESeq() with automatic dispersion fit-type fallback.
# "parametric" requires sufficient spread in gene-wise dispersion estimates and
# fails on very small or very uniform datasets. When it does, we retry with
# "local" (loess-based) and then "mean" (single shared value). If all three
# fail (e.g. all gene-wise dispersions are too tightly clustered), we use the
# gene-wise estimates directly as DESeq2 itself recommends in that case.
deseq_with_fit_fallback <- function(dds, fit_type, size_factor_type, ...) {
  fit_types <- unique(c(fit_type, "local", "mean"))
  last_error <- NULL
  for (ft in fit_types) {
    result <- tryCatch(
      DESeq(dds, fitType = ft, sfType = size_factor_type, ...),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
    if (!is.null(result)) {
      if (ft != fit_type) {
        message(sprintf(
          "DESeq2: fitType='%s' failed (%s); retried with fitType='%s'.",
          fit_type, conditionMessage(last_error), ft
        ))
      }
      return(result)
    }
  }
  # Last resort: use gene-wise dispersion estimates directly (no trend fitting).
  # DESeq2 explicitly recommends this when dispersions are too tightly clustered
  # for any curve to fit, which can happen with small or highly uniform datasets.
  message(sprintf(
    "DESeq2: all fit types (%s) failed (%s); falling back to gene-wise dispersion estimates.",
    paste(fit_types, collapse = ", "),
    conditionMessage(last_error)
  ))
  dds <- estimateSizeFactors(dds, type = size_factor_type)
  dds <- estimateDispersionsGeneEst(dds)
  dispersions(dds) <- mcols(dds)$dispGeneEst
  dds <- nbinomWaldTest(dds)
  dds
}

# Save an MA plot (log2 fold-change vs mean expression) as a PNG.
save_ma_plot <- function(res, path, label) {
  if (!nzchar(trimws(path))) {
    return()
  }
  png(path, width = 800, height = 600, res = 100)
  plotMA(res, main = label)
  dev.off()
}

# === Main script ============================================================
# Guard: only execute the script body when run directly (not when sourced
# by tests).  sys.nframe() == 0L is true only at the top level of a script
# invoked via Rscript; it is FALSE when the file is source()-d.
if (sys.nframe() == 0L) {

# --- Parse command-line arguments -------------------------------------------

counts_path <- get_arg("--counts")
coldata_path <- get_arg("--coldata")
results_path <- get_arg("--results")
size_factors_path <- get_arg("--size-factors")
vst_counts_path <- get_arg("--vst-counts")
summary_path <- get_arg("--summary")
results_names_path <- get_arg("--results-names")
reference_levels_path <- get_arg("--reference-levels")
effect_specs_path <- get_arg("--effect-specs")
fit_type <- get_arg("--fit-type")
size_factor_type <- get_arg("--size-factor-type")
alpha <- as.numeric(get_arg("--alpha"))
cooks_cutoff <- parse_bool(get_arg("--cooks-cutoff"))
independent_filtering <- parse_bool(get_arg("--independent-filtering"))
fixed_effects_formula <- get_arg("--fixed-effects-formula")
test_kind <- tolower(get_arg("--test"))
reduced_formula <- get_arg("--reduced-formula")
ma_plot_path <- get_arg("--ma-plot")

# --- Load DESeq2 and read input data ---------------------------------------

suppressPackageStartupMessages(library("DESeq2"))

counts <- read.table(
  counts_path,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)
coldata <- read.table(
  coldata_path,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

# Ensure integer counts and apply user-specified reference levels.
counts <- round(as.matrix(counts))
storage.mode(counts) <- "integer"
coldata <- apply_reference_levels(coldata, read_list_file(reference_levels_path))

if (!identical(colnames(counts), rownames(coldata))) {
  stop("Sample order mismatch between count matrix and metadata.")
}

# --- Build DESeqDataSet and fit the model -----------------------------------

design_formula <- as.formula(paste("~", fixed_effects_formula))
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = design_formula
)
# Drop features with zero total counts (uninformative, avoids DESeq2 warnings).
dds <- dds[rowSums(counts(dds)) > 0, ]

if (test_kind == "lrt") {
  # Likelihood Ratio Test: compares the full model to a reduced model.
  if (!nzchar(trimws(reduced_formula))) {
    stop("reduced_formula is required when test='lrt'.")
  }
  if (length(read_list_file(effect_specs_path)) > 0) {
    stop("effect_specs are not supported when test='lrt'.")
  }
  dds <- deseq_with_fit_fallback(
    dds,
    fit_type = fit_type,
    size_factor_type = size_factor_type,
    test = "LRT",
    reduced = as.formula(paste("~", reduced_formula))
  )
} else {
  # Default: Wald test for individual coefficients / contrasts.
  dds <- deseq_with_fit_fallback(
    dds,
    fit_type = fit_type,
    size_factor_type = size_factor_type
  )
}

# Write all available coefficient names so the Python side knows what
# can be extracted (e.g. "Intercept", "condition_treated_vs_control").
results_names <- resultsNames(dds)
writeLines(results_names, con = results_names_path)

# --- Extract results for each effect ----------------------------------------

result_frames <- list()
summary_lines <- character()

if (test_kind == "lrt") {
  # LRT produces a single omnibus result (no per-contrast breakdowns).
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
  save_ma_plot(res, ma_plot_path, effect_label)
  result_frames[[effect_id]] <- res_df
  summary_lines <- c(
    summary_lines,
    effect_label,
    capture.output(summary(res)),
    ""
  )
} else {
  # Wald test: extract results for each requested effect spec.
  # If no effect specs are provided, default to all non-intercept coefficients.
  effect_specs <- read_list_file(effect_specs_path)
  if (length(effect_specs) == 0) {
    coefficient_names <- setdiff(results_names, "Intercept")
    effect_specs <- sprintf("coef::%s", coefficient_names)
  }
  if (length(effect_specs) == 0) {
    stop("No non-intercept coefficients were available to extract from the fitted model.")
  }

  # Cache refitted DESeqDataSets for simple effects (same releveling can be
  # shared across multiple numerator levels).
  simple_dds_cache <- list()

  first_effect <- TRUE

  for (spec in effect_specs) {
    if (startsWith(spec, "coef::")) {
      # --- Coefficient extraction: retrieve by name from resultsNames(dds) ---
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
      # --- Pairwise contrast: numerator vs denominator within a factor ------
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
      # --- Simple effect: contrast within a specific level of another factor -
      # e.g. "genotype effect within treatment=compoundA" requires refitting
      # with releveled factors so the interaction contrast is interpretable.
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

      # Refit the model with releveled factors if not already cached.
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
        simple_dds <- deseq_with_fit_fallback(simple_dds, fit_type = fit_type, size_factor_type = size_factor_type)
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

    # Save an MA plot for the first effect only.
    if (first_effect) {
      save_ma_plot(res, ma_plot_path, effect_label)
      first_effect <- FALSE
    }

    result_frames[[spec]] <- res_df
    summary_lines <- c(
      summary_lines,
      res_df$effect_label[[1]],
      capture.output(summary(res)),
      ""
    )
  }
}

# --- Write output files -----------------------------------------------------

# Combine all per-effect result frames into a single table.
res_df <- do.call(rbind, result_frames)

write.table(
  res_df,
  file = results_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Size factors (used by the Python side to compute normalized counts).
size_factors_df <- data.frame(
  sample_id = colnames(dds),
  size_factor = as.numeric(sizeFactors(dds)),
  row.names = NULL,
  check.names = FALSE
)
write.table(
  size_factors_df,
  file = size_factors_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Variance-stabilized counts (used for PCA, distance matrix, and heatmap).
# vst() is faster but requires enough features with mean count > 5 to fit
# its trend; fall back to varianceStabilizingTransformation() for sparse
# datasets.  If both fail (e.g. all dispersions are tightly clustered so no
# trend can be fit), use log2(normalized_counts + 1) as a last resort.
vsd <- tryCatch(
  vst(dds, blind = FALSE, nsub = min(1000, nrow(dds))),
  error = function(e) tryCatch(
    varianceStabilizingTransformation(dds, blind = FALSE),
    error = function(e2) {
      message("VST failed; using log2(normalized counts + 1) for QC visualizations.")
      nc <- counts(dds, normalized = TRUE)
      log2nc <- log2(nc + 1)
      SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = log2nc),
        colData = colData(dds)
      )
    }
  )
)
vst_counts <- as.data.frame(assay(vsd))
vst_counts$feature_id <- rownames(vst_counts)
vst_counts <- vst_counts[
  ,
  c("feature_id", setdiff(colnames(vst_counts), "feature_id"))
]
write.table(
  vst_counts,
  file = vst_counts_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Human-readable summary of each effect (printed by DESeq2's summary()).
if (length(summary_lines) == 0) {
  summary_lines <- "No effects were generated."
}
writeLines(summary_lines, con = summary_path)

} # end if (sys.nframe() == 0L)
