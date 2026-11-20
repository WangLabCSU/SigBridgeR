skip_if_not_installed("scAB")

.scAB_args <- function(phenotype_class = "binary") {
  list(
    phenotype = nb_binary_phenotype(),
    label_type = "scAB",
    phenotype_class = phenotype_class,
    alpha = 0.005,
    alpha_2 = 0.005,
    k_max = 3L,
    cross_k = 2L,
    repeat_times = 2L,
    maxiter = 10L,
    tred = 2L,
    verbose = FALSE,
    seed = 7L
  )
}

test_that("ValidatescABParams resolves defaults", {
  p <- do.call(ValidatescABParams, .scAB_args())

  expect_equal(p$phenotype_class, "binary")
  expect_equal(p$maxiter, 10L)
  expect_equal(p$assay, "RNA")
  expect_false(p$parallel)
  expect_equal(p$seed, 7L)
  expect_true(S7::S7_inherits(p$cache_config, ScreenMethodConfig))
  expect_equal(p$cache_config@method_name, "scAB")
})

test_that("ValidatescABParams supports survival only", {
  p <- do.call(ValidatescABParams, .scAB_args("survival"))
  expect_equal(p$phenotype_class, "survival")

  expect_error(do.call(ValidatescABParams, .scAB_args("continuous")))
})

test_that("ValidatescABParams rejects NA phenotypes and bad inputs", {
  args <- .scAB_args()
  args$phenotype[1L] <- NA
  expect_error(do.call(ValidatescABParams, args))

  args <- .scAB_args()
  args$label_type <- 1L
  expect_error(do.call(ValidatescABParams, args))

  args <- .scAB_args()
  args$maxiter <- "many"
  expect_error(do.call(ValidatescABParams, args))
})

test_that("DoscAB integrates mocked scAB results into the Seurat object", {
  bulk <- nb_bulk(20L, 20L, seed = 11L)
  sc <- nb_seurat(20L, 20L, seed = 12L)
  pheno <- nb_binary_phenotype(20L, seed = 13L)

  local_mocked_bindings(
    create_scAB.v5 = function(...) list(X = matrix(0, 5L, ncol(sc))),
    select_K.optimized = function(...) 2L,
    scAB.optimized = function(...) list(K = 2L),
    findSubset.optimized = function(Object, ...) Object,
    .package = "scAB"
  )

  res <- DoscAB(
    bulk, sc, pheno,
    alpha = 0.005,
    alpha_2 = 0.005,
    k_max = 2L,
    maxiter = 10L,
    verbose = FALSE
  )

  expect_named(res, c("scRNA_data", "scAB_result"))
  expect_true(inherits(res$scRNA_data, "Seurat"))
  expect_equal(ncol(res$scRNA_data), 20L)
  expect_equal(res$scAB_result$K, 2L)
})

test_that("ValidatescABParams forwards cache and parallel settings", {
  args <- c(
    .scAB_args(),
    list(
      parallel = TRUE,
      assay = "SCT",
      load_cache = "/tmp/load",
      save_cache = "/tmp/save"
    )
  )
  p <- do.call(ValidatescABParams, args)

  expect_true(p$parallel)
  expect_equal(p$assay, "SCT")
  expect_equal(p$load_cache, "/tmp/load")
  expect_equal(p$save_cache, "/tmp/save")
})
