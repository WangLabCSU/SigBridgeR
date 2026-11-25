# 11b-SCAnnotate.R ------------------------------------------------------------

# ---- dispatch on unsupported input ------------------------------------------

test_that("SCAnnotate() rejects unsupported input classes", {
  expect_error(
    SCAnnotate(matrix(1:4, nrow = 2L)),
    "Unsupported class",
    class = "rlang_error"
  )
  expect_error(
    SCAnnotate(data.frame(x = 1)),
    "Expected a",
    class = "rlang_error"
  )
})

# ---- SingleR branch ---------------------------------------------------------

test_that("SCAnnotate() dispatches to SingleR with HPCA by default", {
  sc <- nb_seurat(20L, 6L)
  captured <- new.env(parent = emptyenv())

  local_mocked_bindings(
    SingleRAnnotate.Seurat = function(sc, ...) {
      captured$args <- list(...)
      sc
    },
    .package = "SigBridgeR"
  )

  out <- SCAnnotate(sc, method = "SingleR", verbose = FALSE)

  expect_identical(out, sc)
  expect_equal(captured$args$ref, "HPCA")
  expect_false(captured$args$verbose)
})

test_that("SCAnnotate() forwards an explicit reference and partial method names", {
  sc <- nb_seurat(20L, 6L)
  captured <- new.env(parent = emptyenv())

  local_mocked_bindings(
    SingleRAnnotate.Seurat = function(sc, ...) {
      captured$args <- list(...)
      sc
    },
    .package = "SigBridgeR"
  )

  out <- SCAnnotate(
    sc,
    method = "Sing",
    ref = "BlueprintEncode",
    de.n = 10L,
    verbose = FALSE
  )

  expect_identical(out, sc)
  expect_equal(captured$args$ref, "BlueprintEncode")
  expect_equal(captured$args$de.n, 10L)
})

# ---- mLLMCelltype branch ----------------------------------------------------

test_that("SCAnnotate() forwards dots to mLLMCelltype", {
  sc <- nb_seurat(20L, 6L)
  captured <- new.env(parent = emptyenv())

  local_mocked_bindings(
    mLLMCelltypeAnnotate.Seurat = function(sc, ...) {
      captured$args <- list(...)
      sc
    },
    .package = "SigBridgeR"
  )

  out <- SCAnnotate(
    sc,
    method = "mLLMCelltype",
    tissue_name = "Lung Tumor",
    models = "gpt-4o",
    verbose = FALSE
  )

  expect_identical(out, sc)
  expect_equal(captured$args$tissue_name, "Lung Tumor")
  expect_equal(captured$args$models, "gpt-4o")
})

# ---- CellTypist branch ------------------------------------------------------

test_that("SCAnnotate() requires a Python environment for CellTypist", {
  sc <- nb_seurat(20L, 6L)

  local_mocked_bindings(
    ListPyEnv = function(...) data.frame(name = c("envA", "envB")),
    .package = "SigBridgeR"
  )

  expect_error(
    SCAnnotate(sc, method = "CellTypist", verbose = FALSE),
    "python environment",
    class = "rlang_error"
  )
})

test_that("SCAnnotate() forwards CellTypist arguments and defaults", {
  sc <- nb_seurat(20L, 6L)
  captured <- new.env(parent = emptyenv())

  local_mocked_bindings(
    ListPyEnv = function(...) data.frame(name = "fake"),
    CellTypistAnnotate.Seurat = function(sc, ...) {
      captured$args <- list(...)
      sc
    },
    .package = "SigBridgeR"
  )

  out <- SCAnnotate(
    sc,
    method = "CellTypist",
    python = "/fake/python",
    p_thres = 0.8,
    verbose = FALSE
  )

  expect_identical(out, sc)
  expect_true(captured$args$majority_voting)
  expect_equal(captured$args$p_thres, 0.8)

  out <- SCAnnotate(
    sc,
    method = "CellTypist",
    conda = "fake",
    majority_voting = FALSE,
    mode = "prob match",
    verbose = FALSE
  )

  expect_identical(out, sc)
  expect_false(captured$args$majority_voting)
  expect_equal(captured$args$mode, "prob match")
})

# ---- method resolution ------------------------------------------------------

test_that("SCAnnotate() cannot auto-select when several methods remain", {
  sc <- nb_seurat(20L, 6L)

  expect_error(
    SCAnnotate(sc, verbose = FALSE),
    "Cannot auto-find a suitable method",
    class = "rlang_error"
  )
  expect_error(
    SCAnnotate(sc, method = c("SingleR", "CellTypist"), verbose = FALSE),
    "Cannot auto-find a suitable method",
    class = "rlang_error"
  )
})

test_that("SCAnnotate() rejects unknown method names", {
  sc <- nb_seurat(20L, 6L)

  expect_error(
    SCAnnotate(sc, method = "not_a_method", verbose = FALSE),
    "not a valid choice",
    class = "MatchArgError"
  )
})
