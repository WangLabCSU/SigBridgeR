# 50-SingleRAnnotate.R --------------------------------------------------------

skip_if_not_installed("SingleR")
skip_if_not_installed("SingleCellExperiment")

# ---- local fixtures ---------------------------------------------------------

.sr_sce <- function(n_genes = 20L, n_cells = 6L) {
  counts <- matrix(
    seq_len(n_genes * n_cells),
    nrow = n_genes,
    dimnames = list(
      paste0("G", seq_len(n_genes)),
      paste0("C", seq_len(n_cells))
    )
  )
  SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
}

.sr_prediction <- function(test, label = "T cell") {
  data.frame(
    labels = rep(label, ncol(test)),
    delta.next = seq_len(ncol(test)),
    pruned.labels = rep(label, ncol(test)),
    row.names = colnames(test)
  )
}

# ---- input validation -------------------------------------------------------

test_that("SingleRAnnotate() rejects unsupported input classes", {
  expect_error(
    SingleRAnnotate(matrix(1:4, nrow = 2L)),
    "doesn't support",
    class = "rlang_error"
  )
  expect_error(
    SingleRAnnotate(list()),
    "Convert your object to a supported class first",
    class = "rlang_error"
  )
})

test_that("SingleRAnnotate() rejects the 'custom' placeholder reference", {
  expect_error(
    SingleRAnnotate(.sr_sce(), ref = "custom", verbose = FALSE),
    "Please specify the reference dataset",
    class = "rlang_error"
  )
})

test_that("SingleRAnnotate() requires labels for a custom reference", {
  expect_error(
    SingleRAnnotate(.sr_sce(), ref = .sr_sce(), verbose = FALSE),
    "Please specify the `labels`",
    class = "rlang_error"
  )
})

# ---- SingleCellExperiment method -------------------------------------------

test_that("SingleRAnnotate.SingleCellExperiment forwards a custom reference", {
  sce <- .sr_sce()
  ref <- .sr_sce()
  labels <- rep(c("B cell", "T cell"), length.out = ncol(ref))
  captured <- new.env(parent = emptyenv())

  local_mocked_bindings(
    SingleR = function(test, ref, labels, ...) {
      captured$test <- test
      captured$ref <- ref
      captured$labels <- labels
      captured$dots <- list(...)
      .sr_prediction(test)
    },
    .package = "SingleR"
  )

  out <- SingleRAnnotate(
    sce,
    ref = ref,
    labels = labels,
    de.n = 10L,
    verbose = FALSE
  )

  expect_equal(out$labels, rep("T cell", ncol(sce)))
  expect_identical(captured$test, sce)
  expect_identical(captured$ref, ref)
  expect_equal(captured$labels, labels)
  expect_equal(captured$dots$de.n, 10L)
})

test_that("SingleRAnnotate.SingleCellExperiment reports progress when verbose", {
  sce <- .sr_sce()
  ref <- .sr_sce()
  labels <- rep("T cell", ncol(ref))

  local_mocked_bindings(
    SingleR = function(test, ref, labels, ...) .sr_prediction(test),
    .package = "SingleR"
  )

  messages <- capture_messages(
    SingleRAnnotate(sce, ref = ref, labels = labels, verbose = TRUE)
  )
  all_messages <- paste(messages, collapse = "\n")

  expect_match(all_messages, "SingleR", fixed = TRUE)
  expect_match(all_messages, "Annotation finish", fixed = TRUE)
})

test_that("SingleRAnnotate() reports a missing celldex dependency for HPCA", {
  skip_if(
    requireNamespace("celldex", quietly = TRUE),
    "celldex is installed; the missing-dependency path cannot be tested"
  )

  expect_error(
    SingleRAnnotate(.sr_sce(), ref = "HPCA", verbose = FALSE),
    "celldex"
  )
})

test_that("SingleRAnnotate() resolves the HPCA reference through celldex", {
  skip_if_not_installed("celldex")

  sce <- .sr_sce()
  captured <- new.env(parent = emptyenv())
  ref <- .sr_sce()
  ref$label.main <- rep(c("B cell", "T cell"), length.out = ncol(ref))

  local_mocked_bindings(
    HumanPrimaryCellAtlasData = function(ensembl, cell.ont, legacy) {
      captured$celldex <- list(
        ensembl = ensembl,
        cell.ont = cell.ont,
        legacy = legacy
      )
      ref
    },
    .package = "celldex"
  )
  local_mocked_bindings(
    SingleR = function(test, ref, labels, ...) {
      captured$labels <- labels
      .sr_prediction(test)
    },
    .package = "SingleR"
  )

  out <- SingleRAnnotate(
    sce,
    ref = "HPCA",
    ensembl = TRUE,
    cell.ont = "none",
    legacy = TRUE,
    verbose = FALSE
  )

  expect_equal(captured$celldex, list(
    ensembl = TRUE,
    cell.ont = "none",
    legacy = TRUE
  ))
  expect_equal(captured$labels, ref$label.main)
  expect_equal(out$labels, rep("T cell", ncol(sce)))
})

# ---- Seurat method ----------------------------------------------------------

test_that("SingleRAnnotate.Seurat adds the SingleR metadata columns", {
  sc <- nb_seurat(20L, 6L)

  local_mocked_bindings(
    SingleRAnnotate.SingleCellExperiment = function(sc, ...) {
      .sr_prediction(sc, label = "B cell")
    },
    .package = "SigBridgeR"
  )

  out <- suppressWarnings(
    SingleRAnnotate(sc, ref = "HPCA", verbose = FALSE)
  )

  expect_s4_class(out, "Seurat")
  expect_true(all(
    c(
      "SingleR_labels",
      "SingleR_delta_next",
      "SingleR_pruned_labels"
    ) %chin% colnames(out@meta.data)
  ))
  expect_equal(unname(out$SingleR_labels), rep("B cell", ncol(sc)))
  expect_equal(unname(out$SingleR_delta_next), seq_len(ncol(sc)))
  expect_equal(unname(out$SingleR_pruned_labels), rep("B cell", ncol(sc)))
})

# ---- ScreenMethodResult method ---------------------------------------------

test_that("SingleRAnnotate() unwraps a ScreenMethodResult", {
  sc <- nb_seurat(20L, 6L)
  result <- ScreenMethodResult(scRNA_data = sc)
  captured <- new.env(parent = emptyenv())

  local_mocked_bindings(
    SingleRAnnotate.Seurat = function(sc, ...) {
      captured$sc <- sc
      "annotated"
    },
    .package = "SigBridgeR"
  )

  out <- SingleRAnnotate(result)

  expect_equal(out, "annotated")
  expect_identical(captured$sc, sc)
})
