# Unit tests for the helper functions defined in run_deseq2.R.
#
# Run from the repo root:
#   Rscript -e "testthat::test_file('q2_deseq2/r/test_run_deseq2.R')"
#
# Or inside the Docker container (rna-seq env):
#   conda run -n rna-seq Rscript -e \
#     "testthat::test_file('/work/q2-deseq2/q2_deseq2/r/test_run_deseq2.R')"

library(testthat)

# Source the R script so all helper functions are available.  The main body is
# guarded by sys.nframe() == 0L so it will NOT execute when sourced.
script_path <- local({
  # When invoked directly via Rscript the file path is in --file= arg.
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  test_dir <- if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  } else {
    # testthat::test_file() sets cwd to the directory containing the test file.
    getwd()
  }
  normalizePath(file.path(test_dir, "run_deseq2.R"), mustWork = TRUE)
})
source(script_path)

# ============================================================================
# parse_bool
# ============================================================================

test_that("parse_bool recognises truthy strings case-insensitively", {
  expect_true(parse_bool("true"))
  expect_true(parse_bool("TRUE"))
  expect_true(parse_bool("True"))
  expect_true(parse_bool("1"))
  expect_true(parse_bool("yes"))
  expect_true(parse_bool("YES"))
})

test_that("parse_bool treats anything else as FALSE", {
  expect_false(parse_bool("false"))
  expect_false(parse_bool("FALSE"))
  expect_false(parse_bool("0"))
  expect_false(parse_bool("no"))
  expect_false(parse_bool("No"))
  expect_false(parse_bool(""))
  expect_false(parse_bool("maybe"))
})

# ============================================================================
# parse_reference_level
# ============================================================================

test_that("parse_reference_level splits a valid 'column::level' spec", {
  result <- parse_reference_level("condition::control")
  expect_equal(result, c("condition", "control"))
})

test_that("parse_reference_level works with underscores and mixed case", {
  result <- parse_reference_level("sample_type::treated_A")
  expect_equal(result, c("sample_type", "treated_A"))
})

test_that("parse_reference_level errors when separator is absent", {
  expect_error(parse_reference_level("conditioncontrol"))
})

test_that("parse_reference_level errors when a part is empty", {
  expect_error(parse_reference_level("::control"))
  expect_error(parse_reference_level("condition::"))
})

# ============================================================================
# ordered_levels
# ============================================================================

test_that("ordered_levels puts the reference level first", {
  vals <- c("b", "c", "a", "b", "c")
  result <- ordered_levels(vals, "a")
  expect_equal(result[1], "a")
  expect_setequal(result, c("a", "b", "c"))
})

test_that("ordered_levels places reference first even when it sorts last", {
  result <- ordered_levels(c("control", "treated"), "treated")
  expect_equal(result[1], "treated")
  expect_equal(result[2], "control")
})

test_that("ordered_levels returns remaining levels in sorted order", {
  result <- ordered_levels(c("c", "a", "b", "ref"), "ref")
  expect_equal(result, c("ref", "a", "b", "c"))
})

test_that("ordered_levels errors when the reference level is absent", {
  expect_error(ordered_levels(c("a", "b"), "missing"))
})

# ============================================================================
# coerce_coldata
# ============================================================================

test_that("coerce_coldata converts character columns to factors", {
  df <- data.frame(group = c("a", "b", "a"), stringsAsFactors = FALSE)
  result <- coerce_coldata(df)
  expect_true(is.factor(result$group))
})

test_that("coerce_coldata leaves numeric columns unchanged", {
  df <- data.frame(score = c(1.0, 2.5, 3.0), stringsAsFactors = FALSE)
  result <- coerce_coldata(df)
  expect_true(is.numeric(result$score))
})

test_that("coerce_coldata handles a mix of column types", {
  df <- data.frame(
    group = c("a", "b"),
    value = c(10L, 20L),
    stringsAsFactors = FALSE
  )
  result <- coerce_coldata(df)
  expect_true(is.factor(result$group))
  expect_true(is.integer(result$value))
})

# ============================================================================
# apply_reference_levels
# ============================================================================

test_that("apply_reference_levels with no specs just coerces columns to factors", {
  df <- data.frame(group = c("b", "a", "b"), stringsAsFactors = FALSE)
  result <- apply_reference_levels(df, character(0))
  expect_true(is.factor(result$group))
})

test_that("apply_reference_levels sets the requested factor order", {
  df <- data.frame(condition = c("treated", "control", "treated"), stringsAsFactors = FALSE)
  result <- apply_reference_levels(df, "condition::control")
  expect_equal(levels(result$condition)[1], "control")
})

test_that("apply_reference_levels handles multiple specs independently", {
  df <- data.frame(
    condition = c("treated", "control"),
    batch = c("B", "A"),
    stringsAsFactors = FALSE
  )
  result <- apply_reference_levels(df, c("condition::control", "batch::A"))
  expect_equal(levels(result$condition)[1], "control")
  expect_equal(levels(result$batch)[1], "A")
})

test_that("apply_reference_levels errors when the column is absent", {
  df <- data.frame(condition = c("a", "b"), stringsAsFactors = FALSE)
  expect_error(apply_reference_levels(df, "missing_col::a"))
})

test_that("apply_reference_levels errors when the column is numeric", {
  df <- data.frame(score = c(1.0, 2.0), stringsAsFactors = FALSE)
  expect_error(apply_reference_levels(df, "score::1"))
})

test_that("apply_reference_levels errors when the reference level is not present", {
  df <- data.frame(condition = c("a", "b"), stringsAsFactors = FALSE)
  expect_error(apply_reference_levels(df, "condition::missing"))
})

# ============================================================================
# relevel_for_simple_effect
# ============================================================================

make_simple_coldata <- function() {
  data.frame(
    genotype  = c("wt", "mut", "wt", "mut"),
    treatment = c("vehicle", "vehicle", "drug", "drug"),
    stringsAsFactors = FALSE
  )
}

test_that("relevel_for_simple_effect correctly relevels both factors", {
  cd <- make_simple_coldata()
  result <- relevel_for_simple_effect(
    cd,
    factor_name  = "genotype",
    denominator  = "wt",
    within_factor = "treatment",
    within_level  = "drug"
  )
  expect_equal(levels(result$genotype)[1],  "wt")
  expect_equal(levels(result$treatment)[1], "drug")
})

test_that("relevel_for_simple_effect does not modify the original data frame", {
  cd <- make_simple_coldata()
  relevel_for_simple_effect(cd, "genotype", "wt", "treatment", "drug")
  expect_true(is.character(cd$genotype))  # original unchanged
})

test_that("relevel_for_simple_effect errors when factor_name equals within_factor", {
  cd <- make_simple_coldata()
  expect_error(
    relevel_for_simple_effect(cd, "genotype", "wt", "genotype", "mut")
  )
})

test_that("relevel_for_simple_effect errors when factor_name is missing from coldata", {
  cd <- make_simple_coldata()
  expect_error(
    relevel_for_simple_effect(cd, "absent", "wt", "treatment", "drug")
  )
})

test_that("relevel_for_simple_effect errors when within_factor is missing from coldata", {
  cd <- make_simple_coldata()
  expect_error(
    relevel_for_simple_effect(cd, "genotype", "wt", "absent", "drug")
  )
})

test_that("relevel_for_simple_effect errors when factor_name column is numeric", {
  cd <- data.frame(score = c(1, 2), treatment = c("a", "b"), stringsAsFactors = FALSE)
  expect_error(
    relevel_for_simple_effect(cd, "score", "1", "treatment", "a")
  )
})

# ============================================================================
# make_result_frame
# ============================================================================

make_mock_res <- function(n = 5) {
  # A plain data.frame stands in for a DESeqResults object because
  # make_result_frame() only calls as.data.frame() on it.
  set.seed(42)
  data.frame(
    baseMean      = runif(n, 10, 1000),
    log2FoldChange = rnorm(n),
    lfcSE         = runif(n, 0.1, 0.5),
    stat          = rnorm(n),
    pvalue        = runif(n, 0, 1),
    padj          = runif(n, 0, 1),
    row.names     = paste0("gene_", seq_len(n))
  )
}

test_that("make_result_frame returns a data.frame", {
  res <- make_mock_res()
  df <- make_result_frame(res, "coef::condition_treated_vs_control",
                          "treated vs control", "contrast", 'contrast=c("condition","treated","control")')
  expect_true(is.data.frame(df))
})

test_that("make_result_frame prepends required metadata columns", {
  res <- make_mock_res()
  df <- make_result_frame(res, "coef::x", "x", "coefficient", 'name="x"')
  expected_first_cols <- c(
    "effect_id", "effect_label", "effect_kind", "effect_expression",
    "comparison", "test_level", "reference_level", "feature_id"
  )
  expect_equal(colnames(df)[seq_along(expected_first_cols)], expected_first_cols)
})

test_that("make_result_frame propagates effect metadata correctly", {
  res <- make_mock_res()
  df <- make_result_frame(
    res,
    effect_id        = "contrast::condition::treated::control",
    effect_label     = "treated vs control",
    effect_kind      = "contrast",
    effect_expression = 'contrast=c("condition","treated","control")',
    comparison       = "treated vs. control",
    test_level       = "treated",
    reference_level  = "control"
  )
  expect_true(all(df$effect_id       == "contrast::condition::treated::control"))
  expect_true(all(df$test_level      == "treated"))
  expect_true(all(df$reference_level == "control"))
  expect_true(all(df$comparison      == "treated vs. control"))
})

test_that("make_result_frame sets feature_id from rownames", {
  res <- make_mock_res(3)
  df <- make_result_frame(res, "coef::x", "x", "coefficient", 'name="x"')
  expect_equal(sort(df$feature_id), c("gene_1", "gene_2", "gene_3"))
})

test_that("make_result_frame sorts rows so smallest padj comes first", {
  res <- make_mock_res(10)
  res$padj <- c(0.9, 0.01, 0.5, 0.001, 0.3, 0.8, 0.05, NA, 0.7, NA)
  res$pvalue <- res$padj / 2
  df <- make_result_frame(res, "coef::x", "x", "coefficient", 'name="x"')
  non_na_padj <- df$padj[!is.na(df$padj)]
  expect_equal(non_na_padj, sort(non_na_padj))
})

test_that("make_result_frame places NA padj rows after non-NA rows", {
  res <- make_mock_res(4)
  res$padj   <- c(NA, 0.01, NA, 0.5)
  res$pvalue <- c(NA, 0.005, NA, 0.25)
  df <- make_result_frame(res, "coef::x", "x", "coefficient", 'name="x"')
  expect_true(all(is.na(tail(df$padj, 2))))
})

test_that("make_result_frame default comparison/test_level/reference_level are empty", {
  res <- make_mock_res(2)
  df <- make_result_frame(res, "coef::x", "x", "coefficient", 'name="x"')
  expect_true(all(df$comparison      == ""))
  expect_true(all(df$test_level      == ""))
  expect_true(all(df$reference_level == ""))
})

# ============================================================================
# deseq_with_fit_fallback  (requires DESeq2)
# ============================================================================

skip_if_not_installed("DESeq2")
library(DESeq2)

make_minimal_dds <- function(n_genes = 20, seed = 1) {
  set.seed(seed)
  counts <- matrix(
    rnbinom(n_genes * 6, mu = 100, size = 5),
    nrow = n_genes,
    dimnames = list(
      paste0("gene", seq_len(n_genes)),
      paste0("s", seq_len(6))
    )
  )
  coldata <- data.frame(
    condition = factor(rep(c("ctrl", "trt"), each = 3)),
    row.names = paste0("s", seq_len(6))
  )
  DESeqDataSetFromMatrix(counts, colData = coldata, design = ~condition)
}

test_that("deseq_with_fit_fallback returns a fitted DESeqDataSet", {
  dds <- make_minimal_dds()
  result <- deseq_with_fit_fallback(dds, fit_type = "parametric", size_factor_type = "ratio")
  expect_s4_class(result, "DESeqDataSet")
  expect_true(!is.null(sizeFactors(result)))
  expect_true(!is.null(dispersions(result)))
})

test_that("deseq_with_fit_fallback result has non-null Wald test stats", {
  dds <- make_minimal_dds()
  result <- deseq_with_fit_fallback(dds, fit_type = "parametric", size_factor_type = "ratio")
  res <- results(result)
  expect_true("log2FoldChange" %in% names(res))
  expect_true("pvalue" %in% names(res))
})

test_that("deseq_with_fit_fallback falls back to gene-wise estimates on uniform dispersions", {
  # A dataset with extremely low variance (near-identical counts across groups)
  # forces all dispersion fits to fail, triggering the gene-wise fallback.
  set.seed(99)
  n <- 8
  counts <- matrix(
    rep(50L, n * 4),
    nrow = n,
    dimnames = list(paste0("g", seq_len(n)), paste0("s", seq_len(4)))
  )
  # Add tiny jitter so DESeq2 can at least estimate something
  counts[1, ] <- c(50L, 51L, 50L, 51L)
  coldata <- data.frame(
    cond = factor(c("a", "a", "b", "b")),
    row.names = paste0("s", seq_len(4))
  )
  dds <- DESeqDataSetFromMatrix(counts, colData = coldata, design = ~cond)
  # Should not throw even if all standard fit types fail
  expect_no_error(
    deseq_with_fit_fallback(dds, fit_type = "parametric", size_factor_type = "ratio")
  )
})

test_that("deseq_with_fit_fallback works with poscounts size factor type", {
  # poscounts is the recommended method for sparse microbiome data
  set.seed(7)
  n_genes <- 15
  counts <- matrix(
    c(rpois(n_genes * 3, lambda = 50), rpois(n_genes * 3, lambda = 5)),
    nrow = n_genes,
    dimnames = list(
      paste0("asv", seq_len(n_genes)),
      paste0("s", seq_len(6))
    )
  )
  counts[1:3, 4:6] <- 0L   # introduce zeros like a microbiome table
  coldata <- data.frame(
    site = factor(rep(c("gut", "skin"), each = 3)),
    row.names = paste0("s", seq_len(6))
  )
  dds <- DESeqDataSetFromMatrix(counts, colData = coldata, design = ~site)
  result <- deseq_with_fit_fallback(dds, fit_type = "parametric", size_factor_type = "poscounts")
  expect_s4_class(result, "DESeqDataSet")
})
