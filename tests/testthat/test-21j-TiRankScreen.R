skip_if_not_installed("rTIRank")

.tirank_args <- function(phenotype_class = "binary") {
  list(
    label_type = "TiRank",
    phenotype_class = phenotype_class,
    tirank_params = list(n_trials = 2L),
    save_path = "./TiRank_res",
    load_cache = NULL,
    verbose = FALSE,
    seed = 7L
  )
}

test_that("ValidateTiRankParams maps phenotype classes to modes", {
  p <- do.call(ValidateTiRankParams, .tirank_args())

  expect_equal(p$mode, "Classification")
  expect_equal(p$phenotype_class, "binary")
  expect_equal(p$assay, "RNA")
  expect_equal(p$seed, 7L)
  expect_false(p$verbose)
  expect_true(S7::S7_inherits(p$cache_config, ScreenMethodConfig))
  expect_equal(p$cache_config@method_name, "TiRank")

  expect_equal(do.call(ValidateTiRankParams, .tirank_args("survival"))$mode, "Cox")
  expect_equal(do.call(ValidateTiRankParams, .tirank_args("continuous"))$mode, "Regression")
  expect_error(do.call(ValidateTiRankParams, .tirank_args("unknown")))
})

test_that("ValidateTiRankParams resolves cache paths from dots", {
  args <- c(
    .tirank_args(),
    list(
      load_cache = "/tmp/deprecated_cache",
      save_cache = "/tmp/new_cache",
      assay = "SCT"
    )
  )
  p <- do.call(ValidateTiRankParams, args)

  expect_equal(p$load_cache, "/tmp/deprecated_cache")
  expect_equal(p$save_cache, "/tmp/new_cache")
  expect_equal(p$assay, "SCT")

  args <- c(.tirank_args(), list(load_cache = "/tmp/positional"))
  expect_equal(do.call(ValidateTiRankParams, args)$load_cache, "/tmp/positional")
})
