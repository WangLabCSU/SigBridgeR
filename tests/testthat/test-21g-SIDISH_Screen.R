skip_if_not_installed("rSIDISH")

.sidish_args <- function(
  phenotype_class = "survival",
  sc_data = nb_seurat()
) {
  list(
    matched_bulk = nb_bulk(),
    sc_data = sc_data,
    phenotype = nb_survival_phenotype(),
    label_type = "SIDISH",
    phenotype_class = phenotype_class,
    sidish_params = list(phase1_epochs = 5L),
    verbose = FALSE,
    seed = 7L
  )
}

test_that("ValidateSIDISHParams resolves defaults and SIDISH parameters", {
  p <- do.call(ValidateSIDISHParams, .sidish_args())

  expect_equal(p$phenotype_class, "survival")
  expect_equal(p$sidish_params$phase1_epochs, 5L)
  expect_equal(p$sidish_params$device, "cuda")
  expect_equal(p$assay, "RNA")
  expect_equal(p$seed, 7L)
  expect_false(p$verbose)
  expect_true(S7::S7_inherits(p$cache_config, ScreenMethodConfig))
  expect_equal(p$cache_config@method_name, "SIDISH")
})

test_that("ValidateSIDISHParams enforces survival and Seurat input", {
  expect_error(do.call(
    ValidateSIDISHParams,
    .sidish_args(phenotype_class = "binary")
  ))
  expect_error(do.call(
    ValidateSIDISHParams,
    .sidish_args(sc_data = nb_sc_matrix())
  ))
})

test_that("ValidateSIDISHParams warns for deprecated env_params", {
  expect_warning(
    p <- do.call(
      ValidateSIDISHParams,
      c(.sidish_args(), list(env_params = list(env.name = "custom")))
    ),
    "deprecated"
  )
  expect_equal(p$env_params, list(env.name = "custom"))
})

test_that("SIDISHParamSet and SIDISHEnvSet apply overrides", {
  params <- SIDISHParamSet(list(phase1_epochs = 3L, device = "cpu"))
  expect_equal(params$phase1_epochs, 3L)
  expect_equal(params$device, "cpu")
  expect_equal(params$phase2_epochs, 500L)

  expect_warning(env <- SIDISHEnvSet(device = "cpu"), "deprecated")
  expect_equal(env$env.name, "r-reticulate-sidish-cpu")
})
