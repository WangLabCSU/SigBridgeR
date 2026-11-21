skip_if_not_installed("LPSGL")

.lpsgl_args <- function(phenotype_class = "binary") {
  list(
    label_type = "LP_SGL",
    phenotype_class = phenotype_class,
    resolution = 0.6,
    alpha = 0.5,
    nfold = 5L,
    dge_analysis = list(run = FALSE),
    verbose = FALSE,
    seed = 7L
  )
}

test_that("ValidateLPSGLParams maps phenotype classes to families", {
  p <- do.call(ValidateLPSGLParams, .lpsgl_args())
  expect_equal(p$family, "binomial")
  expect_equal(p$phenotype_class, "binary")
  expect_equal(p$resolution, 0.6)
  expect_equal(p$assay, "RNA")
  expect_equal(p$seed, 7L)
  expect_true(S7::S7_inherits(p$cache_config, ScreenMethodConfig))
  expect_equal(p$cache_config@method_name, "LP_SGL")

  expect_equal(do.call(ValidateLPSGLParams, .lpsgl_args("survival"))$family, "cox")
  expect_equal(do.call(ValidateLPSGLParams, .lpsgl_args("continuous"))$family, "gaussian")
})

test_that("DoLP_SGL integrates mocked LP-SGL results into the Seurat object", {
  bulk <- nb_bulk(20L, 20L, seed = 11L)
  sc <- nb_seurat(20L, 20L, seed = 12L)
  pheno <- nb_binary_phenotype(20L, seed = 13L)

  local_mocked_bindings(
    run_leiden_clustering = function(...) factor(rep(1L, ncol(sc))),
    label_cell = function(...) {
      list(
        seurat_obj = sc,
        sgl_fit = list(alpha = 0.5),
        cvfit = list(cvm = 0.1)
      )
    },
    .package = "LPSGL"
  )

  res <- DoLP_SGL(
    bulk,
    sc,
    pheno,
    phenotype_class = "binary",
    resolution = 0.6,
    alpha = 0.5,
    nfold = 3L,
    verbose = FALSE
  )

  expect_named(res, c("scRNA_data", "sgl_fit", "cvfit", "dge_res"))
  expect_true(inherits(res$scRNA_data, "Seurat"))
  expect_equal(ncol(res$scRNA_data), 20L)
  expect_equal(res$sgl_fit$alpha, 0.5)
  expect_null(res$dge_res)
})

test_that("ValidateLPSGLParams handles deprecated family", {
  expect_warning(
    p <- do.call(ValidateLPSGLParams, c(.lpsgl_args(), list(family = "cox"))),
    "deprecated"
  )
  expect_equal(p$phenotype_class, "survival")
  expect_equal(p$family, "cox")
})

test_that("ValidateLPSGLParams rejects bad inputs", {
  expect_error(do.call(
    ValidateLPSGLParams,
    within(.lpsgl_args(), phenotype_class <- "unknown")
  ))
  expect_error(do.call(
    ValidateLPSGLParams,
    within(.lpsgl_args(), dge_analysis <- TRUE)
  ))
})
