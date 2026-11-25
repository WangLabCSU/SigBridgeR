# 51-CellTypistAnnotate.R -----------------------------------------------------

# `PyModule` is a package-level environment rather than a namespace, so it needs
# its own binding helper for tests that replace the Python backend.
local_env_binding <- function(env, name, value, .local_envir = parent.frame()) {
  had <- exists(name, envir = env, inherits = FALSE)
  old <- if (had) get(name, envir = env, inherits = FALSE)
  assign(name, value, envir = env)
  withr::defer(
    {
      if (had) {
        assign(name, old, envir = env)
      } else {
        rm(list = name, envir = env)
      }
    },
    envir = .local_envir
  )
}

.fake_celltypist <- function(annotate = NULL, download = NULL) {
  list(
    annotate_celltypist = annotate %||% function(...) NULL,
    download_models = download %||% function(...) NULL
  )
}

# ---- ListCellTypistModels() -------------------------------------------------

test_that("ListCellTypistModels() lists .pkl models from CELLTYPIST_FOLDER", {
  root <- withr::local_tempdir()
  models_dir <- file.path(root, "data", "models")
  dir.create(models_dir, recursive = TRUE)
  file.create(file.path(
    models_dir,
    c("Immune_All_Low.pkl", "Immune_All_High.pkl", "README.txt")
  ))
  withr::local_envvar(CELLTYPIST_FOLDER = root)

  expect_equal(
    ListCellTypistModels(verbose = FALSE),
    c("Immune_All_High", "Immune_All_Low")
  )
  expect_equal(
    ListCellTypistModels(path = models_dir, verbose = FALSE),
    c("Immune_All_High", "Immune_All_Low")
  )
})

test_that("ListCellTypistModels() prints models only when verbose", {
  root <- withr::local_tempdir()
  models_dir <- file.path(root, "data", "models")
  dir.create(models_dir, recursive = TRUE)
  file.create(file.path(models_dir, "Immune_All_Low.pkl"))
  withr::local_envvar(CELLTYPIST_FOLDER = root)

  out <- capture.output(
    ListCellTypistModels(verbose = TRUE),
    type = "message"
  )
  out <- paste(out, collapse = "\n")
  expect_match(out, "Available CellTypist models", fixed = TRUE)
  expect_match(out, "Immune_All_Low", fixed = TRUE)

  vis <- NULL
  capture.output(
    vis <- withVisible(ListCellTypistModels(verbose = TRUE)),
    type = "message"
  )
  expect_false(vis$visible)
})

test_that("ListCellTypistModels() handles an empty model directory", {
  root <- withr::local_tempdir()
  models_dir <- file.path(root, "data", "models")
  dir.create(models_dir, recursive = TRUE)

  expect_equal(ListCellTypistModels(path = models_dir, verbose = FALSE), character(0L))

  out <- capture.output(
    ListCellTypistModels(path = models_dir, verbose = TRUE),
    type = "message"
  )
  out <- paste(out, collapse = "\n")
  expect_match(out, "No CellTypist models found", fixed = TRUE)
})

test_that("ListCellTypistModels() validates its inputs", {
  root <- withr::local_tempdir()

  expect_error(
    ListCellTypistModels(path = file.path(root, "missing"), verbose = FALSE),
    "not found",
    class = "rlang_error"
  )
  expect_error(ListCellTypistModels(verbose = "yes"), "logical")
  expect_error(ListCellTypistModels(path = 1L, verbose = FALSE), "string")
})

# ---- DownloadCellTypistModels() --------------------------------------------

test_that("DownloadCellTypistModels() appends .pkl and forwards options", {
  captured <- new.env(parent = emptyenv())
  fake <- .fake_celltypist(
    download = function(model, force_update) {
      captured$model <- model
      captured$force_update <- force_update
      invisible(NULL)
    }
  )
  local_env_binding(PyModule, "celltypist", fake)

  DownloadCellTypistModels(model = "Immune_All_Low", force_update = TRUE)
  expect_equal(captured$model, "Immune_All_Low.pkl")
  expect_true(captured$force_update)

  DownloadCellTypistModels(model = c("A.pkl", "B"))
  expect_equal(captured$model, c("A.pkl", "B.pkl"))

  DownloadCellTypistModels()
  expect_null(captured$model)
})

test_that("DownloadCellTypistModels() validates its inputs", {
  expect_error(DownloadCellTypistModels(force_update = "yes"), "logical")
  expect_error(DownloadCellTypistModels(model = 1L), "character")
})

# ---- CellTypistAnnotate() dispatch -----------------------------------------

test_that("CellTypistAnnotate() rejects unsupported input classes", {
  expect_error(
    CellTypistAnnotate(matrix(1:4, nrow = 2L)),
    "doesn't support",
    class = "rlang_error"
  )
})

test_that("CellTypistAnnotate() unwraps a ScreenMethodResult", {
  sc <- nb_seurat(20L, 6L)
  result <- ScreenMethodResult(scRNA_data = sc)
  captured <- new.env(parent = emptyenv())

  local_mocked_bindings(
    CellTypistAnnotate.Seurat = function(sc, ...) {
      captured$sc <- sc
      "annotated"
    },
    .package = "SigBridgeR"
  )

  expect_equal(CellTypistAnnotate(result), "annotated")
  expect_identical(captured$sc, sc)
})

test_that("CellTypistAnnotate() converts a SingleCellExperiment input", {
  skip_if_not_installed("SingleCellExperiment")
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(1:20, nrow = 5L))
  )
  captured <- new.env(parent = emptyenv())

  local_mocked_bindings(
    as.Seurat = function(x, ...) {
      captured$sce <- x
      "seurat-from-sce"
    },
    .package = "Seurat"
  )
  local_mocked_bindings(
    CellTypistAnnotate.Seurat = function(sc, ...) {
      captured$sc <- sc
      "annotated"
    },
    .package = "SigBridgeR"
  )

  expect_equal(CellTypistAnnotate(sce), "annotated")
  expect_identical(captured$sce, sce)
  expect_equal(captured$sc, "seurat-from-sce")
})

# ---- CellTypistAnnotate.Seurat() -------------------------------------------

test_that("CellTypistAnnotate.Seurat adds predictions and stores matrices", {
  skip_if_not_installed("anndataR")
  sc <- nb_seurat(20L, 6L)
  captured <- new.env(parent = emptyenv())
  decision <- matrix(
    1,
    nrow = ncol(sc),
    ncol = 2L,
    dimnames = list(colnames(sc), c("B cell", "T cell"))
  )
  probability <- decision / 2

  fake <- .fake_celltypist(
    annotate = function(adata, ...) {
      captured$adata <- adata
      captured$args <- list(...)
      list(
        data.frame(
          prediction = rep("B cell", ncol(sc)),
          row.names = colnames(sc)
        ),
        decision,
        probability
      )
    }
  )
  local_env_binding(PyModule, "celltypist", fake)
  local_mocked_bindings(
    as_AnnData = function(...) {
      captured$as_ann_data <- list(...)
      "adata"
    },
    .package = "anndataR"
  )

  out <- CellTypistAnnotate(
    sc,
    model = "Immune_All_Low",
    majority_voting = TRUE,
    verbose = FALSE
  )

  expect_s4_class(out, "Seurat")
  expect_true("celltypist_prediction" %chin% colnames(out@meta.data))
  expect_equal(unname(out$celltypist_prediction), rep("B cell", ncol(sc)))

  expect_named(
    out@misc$celltypist,
    c("decision_matrix", "probability_matrix")
  )
  expect_equal(out@misc$celltypist$decision_matrix, decision)
  expect_equal(out@misc$celltypist$probability_matrix, probability)

  expect_equal(captured$adata, "adata")
  expect_equal(captured$as_ann_data$assay_name, "RNA")
  expect_equal(captured$as_ann_data$x_mapping, "data")
  expect_equal(captured$args$model, "Immune_All_Low.pkl")
  expect_true(captured$args$majority_voting)
  expect_equal(captured$args$mode, "best match")
})

test_that("CellTypistAnnotate.Seurat forwards mode, threshold and model path", {
  skip_if_not_installed("anndataR")
  sc <- nb_seurat(20L, 6L)
  captured <- new.env(parent = emptyenv())

  fake <- .fake_celltypist(
    annotate = function(adata, ...) {
      captured$args <- list(...)
      list(
        data.frame(
          prediction = rep("T cell", ncol(sc)),
          row.names = colnames(sc)
        ),
        matrix(1, nrow = ncol(sc), ncol = 2L),
        matrix(0.5, nrow = ncol(sc), ncol = 2L)
      )
    }
  )
  local_env_binding(PyModule, "celltypist", fake)
  local_mocked_bindings(
    as_AnnData = function(...) "adata",
    .package = "anndataR"
  )

  out <- CellTypistAnnotate(
    sc,
    model = "/models/My_Model.pkl",
    mode = "prob match",
    p_thres = 0.9,
    transpose_input = TRUE,
    min_prop = 0.2,
    verbose = FALSE
  )

  expect_s4_class(out, "Seurat")
  expect_equal(captured$args$model, "/models/My_Model.pkl")
  expect_equal(captured$args$mode, "prob match")
  expect_equal(captured$args$p_thres, 0.9)
  expect_true(captured$args$transpose_input)
  expect_equal(captured$args$min_prop, 0.2)
  expect_false(captured$args$majority_voting)
})

test_that("CellTypistAnnotate.Seurat validates its arguments", {
  skip_if_not_installed("anndataR")
  sc <- nb_seurat(20L, 6L)

  expect_error(CellTypistAnnotate(sc, model = 1L), "character")
  expect_error(
    CellTypistAnnotate(sc, mode = "not a mode"),
    "must be one of"
  )
  expect_error(CellTypistAnnotate(sc, download = "yes"), "logical")
  expect_error(CellTypistAnnotate(sc, force_update = "yes"), "logical")
  expect_error(CellTypistAnnotate(sc, majority_voting = "yes"), "logical")
  expect_error(CellTypistAnnotate(sc, p_thres = "high"), "numeric")
  expect_error(CellTypistAnnotate(sc, transpose_input = "yes"), "logical")
  expect_error(CellTypistAnnotate(sc, use_GPU = "yes"), "logical")
  expect_error(CellTypistAnnotate(sc, min_prop = "high"), "numeric")
  expect_error(CellTypistAnnotate(sc, verbose = "yes"), "logical")
  expect_error(CellTypistAnnotate(sc, gene_file = 1L), "character")
  expect_error(CellTypistAnnotate(sc, cell_file = 1L), "character")
})
