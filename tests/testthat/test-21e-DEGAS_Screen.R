skip_if_not_installed("DEGAS")

.degas_args <- function(tmp_dir, phenotype = nb_binary_phenotype()) {
  list(
    select_fraction = 0.05,
    min_thresh = 0.4,
    matched_bulk = nb_bulk(),
    sc_data = nb_seurat(),
    phenotype = phenotype,
    sc_data.pheno_colname = NULL,
    label_type = "DEGAS",
    phenotype_class = "binary",
    tmp_dir = tmp_dir,
    degas_params = list(
      DEGAS.model_type = "BlankClass",
      DEGAS.architecture = "Standard"
    ),
    normality_test_method = "jarque-bera",
    verbose = FALSE
  )
}

test_that("ValidateDEGASParams resolves defaults and common genes", {
  tmp_dir <- withr::local_tempdir()
  p <- do.call(ValidateDEGASParams, .degas_args(tmp_dir))

  expect_equal(p$phenotype_class, "binary")
  expect_equal(p$normality_test_method, "jarque-bera")
  expect_equal(p$assay, "RNA")
  expect_equal(p$tmp_dir, paste0(tmp_dir, "/"))
  expect_equal(p$save_cache, paste0(tmp_dir, "/"))
  expect_equal(length(p$cm_genes), 100L)
  expect_true(inherits(p$t_sc_mat, "Matrix"))
  expect_equal(dim(p$t_sc_mat), c(100L, 100L))
  expect_equal(p$degas_params$DEGAS.architecture, "Standard")
  expect_equal(p$degas_params$DEGAS.model_type, "BlankClass")
  expect_true(S7::S7_inherits(p$cache_config, ScreenMethodConfig))
})

test_that("ValidateDEGASParams rejects unsupported inputs", {
  tmp_dir <- withr::local_tempdir()

  expect_error(do.call(
    ValidateDEGASParams,
    within(.degas_args(tmp_dir), tmp_dir <- NULL)
  ))
  expect_error(do.call(
    ValidateDEGASParams,
    within(.degas_args(tmp_dir), normality_test_method <- "invalid")
  ))
  expect_error(do.call(
    ValidateDEGASParams,
    within(.degas_args(tmp_dir), degas_params <- 1L)
  ))

  args <- .degas_args(tmp_dir)
  args$matched_bulk <- nb_matrix(10L, 10L, prefix_rows = "X")
  expect_error(do.call(ValidateDEGASParams, args))
})

test_that("ValidateDEGASParams warns for deprecated env_params", {
  tmp_dir <- withr::local_tempdir()
  expect_warning(
    p <- do.call(
      ValidateDEGASParams,
      c(.degas_args(tmp_dir), list(env_params = list(env.name = "custom")))
    ),
    "deprecated"
  )
  expect_equal(p$env_params, list(env.name = "custom"))
})

test_that("DEGASModelDetect infers model type from available phenotypes", {
  expect_equal(DEGASModelDetect(NULL, nb_binary_phenotype()), "BlankClass")
  expect_equal(DEGASModelDetect(NULL, nb_survival_phenotype()), "BlankCox")
  expect_equal(DEGASModelDetect("cell_label", NULL), "ClassBlank")
  expect_equal(DEGASModelDetect("cell_label", nb_binary_phenotype()), "ClassClass")
  expect_error(DEGASModelDetect(NULL, NULL))
})

test_that("DEGASParamSet fills defaults and honours user overrides", {
  params <- DEGASParamSet(list(DEGAS.ff_depth = 4L, DEGAS.bag_depth = 2L))

  expect_equal(params$DEGAS.ff_depth, 4L)
  expect_equal(params$DEGAS.bag_depth, 2L)
  expect_equal(params$DEGAS.train_steps, 2000L)
  expect_equal(length(params$DEGAS.model_type), 5L)
})

test_that("DEGASEnvSet is deprecated but returns a configuration", {
  expect_warning(env <- DEGASEnvSet(list(env.name = "custom")), "deprecated")
  expect_equal(env$env.name, "custom")
  expect_equal(env$env.python_verion, "3.9.15")
})
