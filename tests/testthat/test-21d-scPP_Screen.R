skip_if_not_installed("ScPP")
skip_if_not_installed("scPAS")
skip_if_not_installed("dplyr")

.scPP_args <- function(phenotype_class = "binary") {
  list(
    phenotype = nb_binary_phenotype(),
    label_type = "scPP",
    phenotype_class = phenotype_class,
    ref_group = 0L,
    Log2FC_cutoff = 0.585,
    estimate_cutoff = 0.2,
    probs = c(0.2, NULL),
    verbose = FALSE,
    seed = 7L
  )
}

test_that("ValidatescPPParams resolves defaults and lowercases classes", {
  p <- do.call(ValidatescPPParams, .scPP_args("Binary"))

  expect_equal(p$phenotype_class, "binary")
  expect_equal(p$ref_group, 0L)
  expect_equal(p$assay, "RNA")
  expect_false(p$parallel)
  expect_equal(p$seed, 7L)
  expect_true(S7::S7_inherits(p$cache_config, ScreenMethodConfig))
  expect_equal(p$cache_config@method_name, "scPP")

  expect_equal(do.call(ValidatescPPParams, .scPP_args("SURVIVAL"))$phenotype_class, "survival")
  expect_equal(do.call(ValidatescPPParams, .scPP_args("Continuous"))$phenotype_class, "continuous")
})

test_that("ValidatescPPParams validates ranges and phenotypes", {
  expect_error(do.call(
    ValidatescPPParams,
    within(.scPP_args(), probs <- 0.6)
  ))
  expect_error(do.call(
    ValidatescPPParams,
    within(.scPP_args(), ref_group <- "zero")
  ))
  expect_error(do.call(
    ValidatescPPParams,
    within(.scPP_args(), phenotype_class <- "categorical")
  ))

  args <- .scPP_args()
  args$phenotype[1L] <- NA
  expect_error(do.call(ValidatescPPParams, args))
})

test_that("ValidatescPPParams forwards cache, parallel and assay settings", {
  args <- c(
    .scPP_args(),
    list(
      parallel = TRUE,
      assay = "SCT",
      load_cache = "/tmp/load",
      save_cache = "/tmp/save"
    )
  )
  p <- do.call(ValidatescPPParams, args)

  expect_true(p$parallel)
  expect_equal(p$assay, "SCT")
  expect_equal(p$load_cache, "/tmp/load")
  expect_equal(p$save_cache, "/tmp/save")
})
