# Tests for R/64-BenchGeneScore.R
#
# Only the cheap AddModuleScore backend runs for real (10 genes x 8 cells).
# Other backends and the config/S7 dispatch paths are mocked so the file can be
# run standalone via testthat::test_files().

.benchgs_obj <- function() {
  counts <- Matrix::Matrix(
    matrix(seq_len(80L), nrow = 10L, dimnames = list(paste0("g", 1:10), paste0("c", 1:8))),
    sparse = TRUE
  )
  obj <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts))
  Seurat::NormalizeData(obj, verbose = FALSE)
}

.benchgs_pathway <- function() {
  list(First = paste0("g", 1:3), Second = paste0("g", 4:6))
}

test_that("AddModuleScore returns a gene-set-by-cell matrix", {
  obj <- .benchgs_obj()
  meta_before <- colnames(obj[[]])

  res <- BenchGeneScore(
    obj,
    pathway = .benchgs_pathway(),
    method = "AddModuleScore",
    nbin = 2L,
    ctrl = 1L,
    verbose = FALSE
  )

  expect_true(is.matrix(res))
  expect_identical(dim(res), c(2L, 8L))
  expect_identical(rownames(res), names(.benchgs_pathway()))
  expect_identical(colnames(res), colnames(obj))
  expect_true(all(is.finite(res)))
  expect_identical(colnames(obj[[]]), meta_before)
})

test_that("AddModuleScore drops gene sets with no gene in the assay", {
  obj <- .benchgs_obj()
  pathway <- list(
    Valid = paste0("g", 1:3),
    Absent = c("not_a_gene_1", "not_a_gene_2")
  )

  res <- BenchGeneScore(
    obj,
    pathway = pathway,
    method = "AddModuleScore",
    nbin = 2L,
    ctrl = 1L,
    verbose = FALSE
  )

  expect_identical(rownames(res), "Valid")
})

test_that("AddModuleScore aborts when no gene set survives filtering", {
  obj <- .benchgs_obj()

  expect_error(
    BenchGeneScore(
      obj,
      pathway = list(Absent = c("not_a_gene_1", "not_a_gene_2")),
      method = "AddModuleScore",
      nbin = 2L,
      ctrl = 1L,
      verbose = FALSE
    ),
    "None of the gene sets"
  )
})

test_that("BenchGeneScore.Seurat validates the backend name", {
  obj <- .benchgs_obj()

  expect_error(
    BenchGeneScore(obj, .benchgs_pathway(), method = "bogus"),
    "must be one of"
  )
  expect_error(
    BenchGeneScore(list()),
    "no applicable method for 'BenchGeneScore'"
  )
})

test_that("BenchGeneScore.Seurat dispatches to the selected backend", {
  obj <- .benchgs_obj()
  seen <- new.env(parent = emptyenv())

  local_mocked_bindings(
    check_installed = function(...) NULL,
    .package = "rlang"
  )
  local_mocked_bindings(
    BenchGeneScoreSeurat.ssGSEA = function(obj, ...) {
      seen$backend <- "ssGSEA"
      "ss"
    },
    BenchGeneScoreSeurat.AUCell = function(obj, ...) {
      seen$backend <- "AUCell"
      "au"
    },
    BenchGeneScoreSeurat.UCell = function(obj, ...) {
      seen$backend <- "UCell"
      "uc"
    },
    BenchGeneScoreSeurat.AddModuleScore = function(obj, ...) {
      seen$backend <- "AddModuleScore"
      "am"
    },
    .package = "SigBridgeR"
  )

  expect_identical(
    vapply(
      c("ssGSEA", "AUCell", "UCell", "AddModuleScore"),
      function(method) BenchGeneScore(obj, .benchgs_pathway(), method = method),
      character(1L)
    ),
    c(ssGSEA = "ss", AUCell = "au", UCell = "uc", AddModuleScore = "am")
  )
  expect_identical(seen$backend, "AddModuleScore")
})

test_that("BenchGeneScore.ScreenMethodResult unwraps scRNA_data", {
  obj <- .benchgs_obj()
  result <- ScreenMethodResult(scRNA_data = obj)
  seen <- new.env(parent = emptyenv())

  local_mocked_bindings(
    BenchGeneScore.Seurat = function(obj, ...) {
      seen$obj <- obj
      invisible(NULL)
    },
    .package = "SigBridgeR"
  )

  BenchGeneScore(result, pathway = .benchgs_pathway(), method = "AddModuleScore")

  expect_identical(seen$obj, obj)
})

test_that("BenchGeneScore.BenchGeneScoreConfig screens methods before scoring", {
  obj <- .benchgs_obj()
  config <- BenchGeneScoreConfig(
    matched_bulk = matrix(seq_len(4L), nrow = 2L, dimnames = list(c("g1", "g2"), c("s1", "s2"))),
    sc_data = obj,
    phenotype = stats::setNames(rep(c(0, 1), 4L), paste0("c", 1:8)),
    label_type = "{method_name}",
    phenotype_class = "binary",
    screen_method = "scPP",
    method_param = MethodParamSet(),
    bench_config = list(method = "AddModuleScore", nbin = 2L, ctrl = 1L)
  )
  seen <- new.env(parent = emptyenv())
  seen$n <- 0L

  local_mocked_bindings(
    Screen = function(...) {
      seen$n <- seen$n + 1L
      ScreenMethodResult(scRNA_data = obj)
    },
    BenchGeneScore.Seurat = function(obj, ...) {
      seen$obj <- obj
      invisible(NULL)
    },
    .package = "SigBridgeR"
  )

  BenchGeneScore(config, verbose = FALSE)

  expect_equal(seen$n, 1L)
  expect_identical(seen$obj, obj)
})
