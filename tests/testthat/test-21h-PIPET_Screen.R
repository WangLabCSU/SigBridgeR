skip_if_not_installed("PIPET")

.pipet_args <- function(phenotype_class = "binary") {
  list(
    phenotype_class = phenotype_class,
    only_pos_marker = FALSE,
    group = NULL,
    discretize_method = "kmeans",
    cutoff = NULL,
    label_type = "PIPET",
    marker_finder = "limma",
    log2FC = 1,
    p_adjust = 0.05,
    show_log2FC = TRUE,
    freq_counts = NULL,
    normalize = TRUE,
    scale = TRUE,
    nPerm = 10L,
    distance = "cosine",
    verbose = FALSE,
    seed = 7L
  )
}

test_that("ValidatePIPETParams resolves defaults", {
  p <- do.call(ValidatePIPETParams, .pipet_args())

  expect_equal(p$phenotype_class, "binary")
  expect_equal(p$marker_finder, "limma")
  expect_equal(p$distance, "cosine")
  expect_equal(p$assay, "RNA")
  expect_equal(p$seed, 7L)
  expect_false(p$parallel)
  expect_true(S7::S7_inherits(p$cache_config, ScreenMethodConfig))
  expect_equal(p$cache_config@method_name, "PIPET")
})

test_that("ValidatePIPETParams validates categorical choices", {
  expect_equal(
    do.call(ValidatePIPETParams, .pipet_args("survival"))$phenotype_class,
    "survival"
  )
  expect_error(do.call(
    ValidatePIPETParams,
    within(.pipet_args(), marker_finder <- "unknown")
  ))
  expect_error(do.call(
    ValidatePIPETParams,
    within(.pipet_args(), distance <- "unknown")
  ))
  expect_error(do.call(
    ValidatePIPETParams,
    within(.pipet_args(), phenotype_class <- "unknown")
  ))
  expect_error(do.call(
    ValidatePIPETParams,
    within(.pipet_args(), nPerm <- 10)
  ))
  expect_error(do.call(
    ValidatePIPETParams,
    within(.pipet_args(), group <- 1L)
  ))
})

test_that("DoPIPET attaches mocked marker-based labels to the Seurat object", {
  bulk <- nb_bulk(20L, 20L, seed = 11L)
  sc <- nb_seurat(20L, 20L, seed = 12L)
  pheno <- nb_binary_phenotype(20L, seed = 13L)

  local_mocked_bindings(
    AdaptPheno = function(...) data.frame(class = rep(c(0L, 1L), 10L)),
    Create_Markers2 = function(...) {
      data.frame(gene = c("G1", "G2"), log2FoldChange = c(2, -2))
    },
    PIPET = function(...) {
      data.frame(
        PIPET_label = rep("Positive", ncol(sc)),
        row.names = colnames(sc)
      )
    },
    .package = "PIPET"
  )

  res <- DoPIPET(
    bulk, sc, pheno,
    phenotype_class = "binary",
    nPerm = 5L,
    marker_finder = "limma",
    verbose = FALSE
  )

  expect_named(res, c("scRNA_data", "markers"))
  expect_true(inherits(res$scRNA_data, "Seurat"))
  expect_equal(ncol(res$scRNA_data), 20L)
  expect_equal(unname(res$scRNA_data$PIPET_label), rep("Positive", 20L))
  expect_equal(nrow(res$markers), 2L)
})

test_that("ValidatePIPETParams forwards cache and DESeq2 marker finder", {
  args <- .pipet_args("continuous")
  args$marker_finder <- "DESeq2"
  args$distance <- "spearman"
  args$parallel <- TRUE
  args$assay <- "SCT"
  args$load_cache <- "/tmp/load"
  args$save_cache <- "/tmp/save"
  p <- do.call(ValidatePIPETParams, args)

  expect_equal(p$phenotype_class, "continuous")
  expect_equal(p$marker_finder, "DESeq2")
  expect_equal(p$distance, "spearman")
  expect_true(p$parallel)
  expect_equal(p$assay, "SCT")
  expect_equal(p$load_cache, "/tmp/load")
  expect_equal(p$save_cache, "/tmp/save")
})
