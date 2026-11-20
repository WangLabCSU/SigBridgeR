# Tests for R/62-BenchGSEA.R
#
# The FGSEA path is exercised on a 10-gene x 12-cell object with 10
# permutations; dispatch and input validation are mocked so the file stays
# cheap enough to run standalone via testthat::test_files().

.benchgsea_obj <- function() {
  counts <- Matrix::Matrix(
    matrix(seq_len(120L), nrow = 10L,
           dimnames = list(paste0("g", 1:10), paste0("c", 1:12))),
    sparse = TRUE
  )
  obj <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts))
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj$truth <- rep(c("Positive", "Negative"), 6L)
  obj$Scissor <- rep(c("Positive", "Negative", "Negative", "Negative"), 3L)
  obj$scPAS <- rep(c("Positive", "Positive", "Negative", "Negative"), 3L)
  obj
}

.benchgsea_pathway <- function() {
  list(first = paste0("g", 1:3), second = paste0("g", 4:6))
}

test_that("BenchGSEA.default reports unsupported classes", {
  expect_error(BenchGSEA(list()), "not implemented for class")
})

test_that("BenchGSEA rejects reserved FindMarkers arguments", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")
  obj <- .benchgsea_obj()

  expect_error(
    BenchGSEA(
      obj,
      group = obj$truth,
      pathway = .benchgsea_pathway(),
      ident.1 = "Positive"
    ),
    "reserved for BenchGSEA"
  )
  expect_error(
    BenchGSEA(
      obj,
      group = obj$truth,
      pathway = .benchgsea_pathway(),
      group.by = "truth"
    ),
    "reserved for BenchGSEA"
  )
})

test_that("BenchGSEA aborts when no screening method column is present", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")
  obj <- .benchgsea_obj()
  obj[["Scissor"]] <- NULL
  obj[["scPAS"]] <- NULL

  expect_error(
    BenchGSEA(obj, group = obj$truth, pathway = .benchgsea_pathway()),
    "has no screen method results"
  )
})

test_that("BenchGSEA.ScreenMethodResult dispatches to the Seurat method", {
  obj <- .benchgsea_obj()
  result <- ScreenMethodResult(scRNA_data = obj)
  seen <- new.env(parent = emptyenv())

  local_mocked_bindings(
    BenchGSEA.Seurat = function(obj, ...) {
      seen$obj <- obj
      invisible(NULL)
    },
    .package = "SigBridgeR"
  )

  BenchGSEA(result, group = obj$truth, pathway = .benchgsea_pathway())

  expect_identical(seen$obj, result)
})

test_that("BenchGSEA.BenchGSEAConfig screens each method before GSEA", {
  obj <- .benchgsea_obj()
  config <- BenchGSEAConfig(
    matched_bulk = matrix(seq_len(4L), nrow = 2L,
                          dimnames = list(c("g1", "g2"), c("s1", "s2"))),
    sc_data = obj,
    phenotype = rep(c(0, 1), 6L) |> stats::setNames(paste0("c", 1:12)),
    label_type = "{method_name}",
    phenotype_class = "binary",
    screen_method = "scPP",
    method_param = MethodParamSet(),
    bench_config = list(group = obj$truth, pathway = .benchgsea_pathway())
  )
  seen <- new.env(parent = emptyenv())
  seen$n <- 0L

  local_mocked_bindings(
    Screen = function(...) {
      seen$n <- seen$n + 1L
      ScreenMethodResult(scRNA_data = obj)
    },
    BenchGSEA.Seurat = function(obj, ...) {
      seen$obj <- obj
      invisible(NULL)
    },
    .package = "SigBridgeR"
  )

  BenchGSEA(config, verbose = FALSE)

  expect_equal(seen$n, 1L)
  expect_identical(seen$obj, obj)
})

test_that("BenchGSEA binarizes screening labels before calling FindMarkers", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")
  obj <- .benchgsea_obj()
  seen <- new.env(parent = emptyenv())

  local_mocked_bindings(
    FindMarkers = function(object, ident.1, ident.2, group.by, ...) {
      seen$meta <- object[[]]
      seen$ident <- c(ident.1, ident.2)
      seen$group_by <- group.by
      stop("stop after capturing arguments")
    },
    .package = "Seurat"
  )

  expect_error(
    BenchGSEA(obj, group = obj$truth, pathway = .benchgsea_pathway(), verbose = FALSE),
    "stop after capturing arguments"
  )
  expect_identical(seen$ident, c("Positive", "Other"))
  expect_identical(seen$group_by, "Scissor")
  expect_setequal(unique(seen$meta$Scissor), c("Positive", "Other"))
  expect_setequal(unique(seen$meta$scPAS), c("Positive", "Other"))
  expect_setequal(unique(seen$meta$truth), c("Positive", "Negative"))
})

test_that("BenchGSEA forwards extra arguments to FindMarkers", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")
  obj <- .benchgsea_obj()
  seen <- new.env(parent = emptyenv())

  local_mocked_bindings(
    FindMarkers = function(object, ident.1, ident.2, group.by, ...) {
      seen$dots <- list(...)
      stop("stop after capturing arguments")
    },
    .package = "Seurat"
  )

  expect_error(
    BenchGSEA(
      obj,
      group = obj$truth,
      pathway = .benchgsea_pathway(),
      p_value_threshold = 0.05,
      logfc.threshold = 0,
      verbose = FALSE
    ),
    "stop after capturing arguments"
  )
  expect_equal(seen$dots$logfc.threshold, 0)
})
