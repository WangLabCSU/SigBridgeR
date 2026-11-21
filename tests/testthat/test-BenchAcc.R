# Tests for R/61-BenchAcc.R
#
# All fixtures are tiny (4 genes x 6 cells) and the metric computation is the
# package's own cheap C++ routine, so this file is safe to run on its own with
# testthat::test_files().

.benchacc_obj <- function() {
  counts <- Matrix::Matrix(
    matrix(seq_len(24L), nrow = 4L, dimnames = list(paste0("g", 1:4), paste0("c", 1:6))),
    sparse = TRUE
  )
  obj <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts))
  obj$truth <- c("Positive", "Negative", "Positive", "Negative", "Positive", "Negative")
  obj$Scissor <- c("Positive", "Negative", "Positive", "Negative", "Negative", "Negative")
  obj$scPAS <- c("Positive", "Positive", "Negative", "Negative", "Negative", "Negative")
  obj
}

.benchacc_config <- function(obj, method = "scPP") {
  BenchAccConfig(
    matched_bulk = matrix(seq_len(4L), nrow = 2L, dimnames = list(c("g1", "g2"), c("s1", "s2"))),
    sc_data = obj,
    phenotype = c(c1 = 0, c2 = 1, c3 = 1, c4 = 0, c5 = 0, c6 = 1),
    label_type = "{method_name}",
    phenotype_class = "binary",
    screen_method = method,
    method_param = MethodParamSet(),
    bench_config = obj$truth
  )
}

test_that("BenchAcc computes all eight binary metrics per screening method", {
  obj <- .benchacc_obj()

  res <- BenchAcc(obj, bench_config = "truth", verbose = FALSE, nthreads = 1L)

  expect_type(res, "list")
  expect_length(res, 2L)
  expect_setequal(names(res), c("Scissor", "scPAS"))
  expect_identical(dim(res$Scissor), c(1L, 8L))
  expect_identical(
    colnames(res$Scissor),
    c("TPR", "FPR", "TNR", "FNR", "F1", "Accuracy", "Precision", "Recall")
  )
  expect_identical(rownames(res$Scissor), "Scissor")
  expect_equal(
    res$Scissor[1L, ],
    c(TPR = 2 / 3, FPR = 0, TNR = 1, FNR = 1 / 3, F1 = 0.8, Accuracy = 5 / 6, Precision = 1, Recall = 2 / 3)
  )
  expect_equal(
    res$scPAS[1L, ],
    c(TPR = 1 / 3, FPR = 1 / 3, TNR = 2 / 3, FNR = 2 / 3, F1 = 0.4, Accuracy = 0.5, Precision = 1 / 2, Recall = 1 / 3)
  )
})

test_that("BenchAcc accepts a per-cell ground-truth vector", {
  obj <- .benchacc_obj()

  from_column <- BenchAcc(obj, "truth", verbose = FALSE)
  from_vector <- BenchAcc(obj, obj$truth, verbose = FALSE)

  expect_equal(from_vector, from_column)
})

test_that("BenchAcc detects screening method columns case-insensitively", {
  obj <- .benchacc_obj()
  obj[["Scissor"]] <- NULL
  obj[["scPAS"]] <- NULL
  obj$SCISSOR <- .benchacc_obj()$Scissor
  obj$scpas <- .benchacc_obj()$scPAS

  res <- BenchAcc(obj, bench_config = "truth", verbose = FALSE)

  expect_setequal(names(res), c("SCISSOR", "scpas"))
  expect_identical(rownames(res$SCISSOR), "SCISSOR")
})

test_that("BenchAcc reports progress when verbose", {
  obj <- .benchacc_obj()

  expect_message(
    BenchAcc(obj, bench_config = "truth", verbose = TRUE),
    "Calculating metrics"
  )
})

test_that("BenchAcc validates the ground-truth specification", {
  obj <- .benchacc_obj()
  obj$all_negative <- rep("Negative", 6L)

  expect_error(BenchAcc(obj, character(), verbose = FALSE), "not specified")
  expect_error(
    BenchAcc(obj, "missing_truth", verbose = FALSE),
    "is not a valid column name"
  )
  expect_error(
    BenchAcc(obj, "all_negative", verbose = FALSE),
    "doesn't have a `Positive` label"
  )
  expect_error(
    BenchAcc(obj, rep("Negative", 6L), verbose = FALSE),
    "doesn't have a `Positive` label"
  )
  expect_error(
    BenchAcc(obj, c("Positive", "Negative"), verbose = FALSE),
    "doesn't have the same length as cells"
  )
})

test_that("BenchAcc aborts when no screening method column is present", {
  obj <- .benchacc_obj()
  obj[["Scissor"]] <- NULL
  obj[["scPAS"]] <- NULL

  expect_error(
    BenchAcc(obj, "truth", verbose = FALSE),
    "has no screen method results"
  )
})

test_that("BenchAcc rejects invalid labels in predicted columns", {
  obj <- .benchacc_obj()
  obj$Scissor[1L] <- "Maybe"

  expect_error(
    BenchAcc(obj, "truth", verbose = FALSE),
    "contains invalid value or NA"
  )
})

test_that("BenchAcc.ScreenMethodResult unwraps scRNA_data", {
  obj <- .benchacc_obj()
  result <- ScreenMethodResult(scRNA_data = obj)

  res <- BenchAcc(result, bench_config = "truth", verbose = FALSE)

  expect_equal(res, BenchAcc(obj, "truth", verbose = FALSE))
})

test_that("BenchAcc.BenchAccConfig screens every configured method once", {
  obj <- .benchacc_obj()
  config <- .benchacc_config(obj)
  seen <- new.env(parent = emptyenv())
  seen$n <- 0L

  local_mocked_bindings(
    Screen = function(...) {
      seen$n <- seen$n + 1L
      ScreenMethodResult(scRNA_data = obj)
    },
    .package = "SigBridgeR"
  )

  res <- BenchAcc(config, verbose = FALSE, nthreads = 1L)

  expect_equal(seen$n, 1L)
  expect_setequal(names(res), c("Scissor", "scPAS"))
  expect_equal(res$Scissor[1L, "Accuracy"], 5 / 6)
})

test_that("BenchAcc.default reports unsupported classes", {
  expect_error(BenchAcc(list()), "not implemented for class")
})
