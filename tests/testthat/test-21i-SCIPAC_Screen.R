skip_if_not_installed("SCIPAC")
skip_if_not_installed("dplyr")

.scipac_args <- function(
  phenotype_class = "binary",
  sc_data = nb_seurat()
) {
  list(
    label_type = "SCIPAC",
    phenotype_class = phenotype_class,
    hvg = 50L,
    do_pca_sc = FALSE,
    n_pc = 10L,
    sc_batch_col = NULL,
    resolution = 1,
    ela_net_alpha = 0.4,
    bt_size = 5L,
    ncore = 1L,
    ci_alpha = 0.05,
    nfold = 2L,
    sc_data = sc_data,
    verbose = FALSE,
    seed = 7L
  )
}

test_that("ValidateSCIPACParams resolves defaults and families", {
  p <- do.call(ValidateSCIPACParams, .scipac_args())

  expect_equal(p$family, "binomial")
  expect_equal(p$phenotype_class, "binary")
  expect_equal(p$hvg, 50L)
  expect_equal(p$assay, "RNA")
  expect_equal(p$seed, 7L)
  expect_false(p$verbose)
  expect_true(S7::S7_inherits(p$cache_config, ScreenMethodConfig))
  expect_equal(p$cache_config@method_name, "SCIPAC")

  expect_equal(do.call(ValidateSCIPACParams, .scipac_args("survival"))$family, "cox")
  expect_equal(do.call(ValidateSCIPACParams, .scipac_args("continuous"))$family, "gaussian")
})

test_that("ValidateSCIPACParams rejects bad inputs", {
  expect_error(do.call(
    ValidateSCIPACParams,
    within(.scipac_args(), phenotype_class <- "unknown")
  ))
  expect_error(do.call(
    ValidateSCIPACParams,
    .scipac_args(sc_data = nb_sc_matrix())
  ))
  expect_error(do.call(
    ValidateSCIPACParams,
    within(.scipac_args(), hvg <- 50)
  ))
})

test_that("OverLapSCBulk aligns bulk and single-cell genes", {
  # SCIPAC::obtain.preprocessed.data() hard-codes min.features = 200, so the
  # overlapping gene set must exceed that threshold for cells to survive.
  bulk <- nb_bulk(300L, 50L)
  sc <- nb_seurat(300L, 50L)

  out <- OverLapSCBulk(sc_data = sc, bulk = bulk, hvg = 20L)

  expect_named(out, c("sc_data_preprocessed", "bulk_preprocessed"))
  expect_equal(nrow(out$bulk_preprocessed), 20L)
  expect_equal(ncol(out$bulk_preprocessed), 50L)

  other_bulk <- nb_matrix(10L, 10L, prefix_rows = "X")
  expect_error(OverLapSCBulk(sc_data = sc, bulk = other_bulk, hvg = 5L))
})

# DoSCIPAC - wrapper behaviour ---------------------------------------------

test_that("DoSCIPAC aborts when sc_batch_col matches no metadata column", {
  # Same fixtures as the OverLapSCBulk test: the abort sits after gene
  # alignment, so the shared gene set must clear min.features = 200.
  bulk <- nb_bulk(300L, 50L)
  sc <- nb_seurat(300L, 50L)
  sc$batch_id <- rep(c("Batch1", "Batch2"), length.out = ncol(sc))

  pheno <- nb_binary_phenotype(50L)

  expect_error(
    DoSCIPAC(
      matched_bulk = bulk,
      sc_data = sc,
      phenotype = pheno,
      phenotype_class = "binary",
      hvg = 20L,
      sc_batch_col = "batch_typo",
      verbose = FALSE,
      seed = 7L
    ),
    "no batch information found"
  )
})
