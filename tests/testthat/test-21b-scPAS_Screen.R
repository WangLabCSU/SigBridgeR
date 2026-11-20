skip_if_not_installed("scPAS")
skip_if_not_installed("dplyr")

# ValidatescPASParams - defaults & family mapping --------------------------

.scPAS_args <- function(label_type = "scPAS", phenotype_class = "binary") {
  list(
    phenotype = nb_binary_phenotype(),
    label_type = label_type,
    phenotype_class = phenotype_class,
    assay = "RNA",
    imputation = FALSE,
    imputation_method = "KNN",
    nfeature = 50L,
    alpha = 0.01,
    cutoff = 0.2,
    network_class = "SC",
    permutation_times = 20L,
    FDR_threshold = 0.05,
    independent = TRUE,
    verbose = FALSE,
    seed = 7L
  )
}

test_that("ValidatescPASParams resolves binary defaults", {
  p <- do.call(ValidatescPASParams, .scPAS_args())

  expect_equal(p$family, "binomial")
  expect_equal(as.character(p$label_type), c("scPAS_Positive", "scPAS_Negative"))
  expect_equal(p$network_class, "SC")
  expect_equal(p$imputation_method, "KNN")
  expect_equal(p$permutation_times, 20L)
  expect_equal(p$seed, 7L)
  expect_false(p$verbose)
  expect_true(S7::S7_inherits(p$cache_config, ScreenMethodConfig))
  expect_equal(p$cache_config@method_name, "scPAS")
})

test_that("ValidatescPASParams maps survival and continuous phenotypes", {
  p_surv <- do.call(ValidatescPASParams, .scPAS_args(phenotype_class = "survival"))
  expect_equal(p_surv$family, "cox")

  p_cont <- do.call(ValidatescPASParams, .scPAS_args(phenotype_class = "continuous"))
  expect_equal(p_cont$family, "gaussian")
})

test_that("ValidatescPASParams validates choices and numeric ranges", {
  expect_error(do.call(
    ValidatescPASParams,
    within(.scPAS_args(), network_class <- "invalid")
  ))
  expect_error(do.call(
    ValidatescPASParams,
    within(.scPAS_args(), alpha <- 1.5)
  ))
  expect_error(do.call(
    ValidatescPASParams,
    within(.scPAS_args(), {
      imputation <- TRUE
      imputation_method <- "invalid"
    })
  ))
  expect_error(do.call(
    ValidatescPASParams,
    within(.scPAS_args(), label_type <- 1L)
  ))
})

test_that("ValidatescPASParams accepts bulk network and ALRA imputation", {
  args <- .scPAS_args()
  args$network_class <- "bulk"
  args$imputation <- TRUE
  args$imputation_method <- "ALRA"
  p <- do.call(ValidatescPASParams, args)

  expect_equal(p$network_class, "bulk")
  expect_true(p$imputation)
  expect_equal(p$imputation_method, "ALRA")
})

test_that("ValidatescPASParams warns on deprecated family", {
  expect_warning(
    p <- do.call(ValidatescPASParams, c(.scPAS_args(), list(family = "gaussian"))),
    "deprecated"
  )
  expect_equal(p$family, "gaussian")
  expect_equal(p$phenotype_class, "continuous")
})
