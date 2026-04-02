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

save_ma_plot <- function(res, path, label) {
  if (!nzchar(trimws(path))) {
    return()
  }
  png(path, width = 800, height = 600, res = 100)
  plotMA(res, main = label)
  dev.off()
}

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
alpha <- as.numeric(get_arg("--alpha"))
cooks_cutoff <- parse_bool(get_arg("--cooks-cutoff"))
independent_filtering <- parse_bool(get_arg("--independent-filtering"))
fixed_effects_formula <- get_arg("--fixed-effects-formula")
test_kind <- tolower(get_arg("--test"))
reduced_formula <- get_arg("--reduced-formula")
ma_plot_path <- get_arg("--ma-plot")

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
    fitType = fit_type,
    sfType = "poscounts"
  )
} else {
  dds <- DESeq(dds, fitType = fit_type, sfType = "poscounts")
}

results_names <- resultsNames(dds)
writeLines(results_names, con = results_names_path)

result_frames <- list()
summary_lines <- character()

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
  save_ma_plot(res, ma_plot_path, effect_label)
  result_frames[[effect_id]] <- res_df
  summary_lines <- c(
    summary_lines,
    effect_label,
    capture.output(summary(res)),
    ""
  )
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

  first_effect <- TRUE

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
        simple_dds <- DESeq(simple_dds, fitType = fit_type, sfType = "poscounts")
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

res_df <- do.call(rbind, result_frames)

write.table(
  res_df,
  file = results_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

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

vsd <- vst(dds, blind = FALSE, nsub = min(1000, nrow(dds)))
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

if (length(summary_lines) == 0) {
  summary_lines <- "No effects were generated."
}
writeLines(summary_lines, con = summary_path)
