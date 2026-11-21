skip_if_not_installed("Scissor")

# ValidateScissorParams - defaults & family mapping ------------------------

test_that("ValidateScissorParams resolves binary defaults", {
  p <- ValidateScissorParams(
    phenotype = nb_binary_phenotype(),
    label_type = "scissor",
    phenotype_class = "binary",
    alpha = 0.05,
    cutoff = 0.2,
    reliability_test = list(),
    cell_evaluation = list(),
    verbose = FALSE
  )

  expect_equal(p$family, "binomial")
  expect_equal(p$phenotype_class, "binary")
  expect_equal(as.character(p$tag), c("scissor_Negative", "scissor_Positive"))
  expect_equal(p$reliability_test, list(run = FALSE, n = 10L, nfold = 10L))
  expect_equal(
    p$cell_evaluation,
    list(
      run = FALSE,
      benchmark_data = "path_to_file.RData",
      FDR = 0.05,
      bootstrap_n = 100L
    )
  )
  expect_equal(p$assay, "RNA")
  expect_false(p$verbose)
  expect_true(S7::S7_inherits(p$cache_config, ScreenMethodConfig))
  expect_equal(p$cache_config@method_name, "Scissor")
})

test_that("ValidateScissorParams maps survival and continuous phenotypes", {
  pheno <- nb_binary_phenotype()

  p_surv <- ValidateScissorParams(
    phenotype = pheno,
    label_type = "OS",
    phenotype_class = "survival",
    alpha = 0.05,
    cutoff = 0.2,
    reliability_test = list(),
    cell_evaluation = list()
  )
  expect_equal(p_surv$family, "cox")
  expect_equal(as.character(p_surv$tag), c("OS_Negative", "OS_Positive"))

  cont <- c(a = 1, b = 2, c = 1, d = 3)
  p_cont <- ValidateScissorParams(
    phenotype = cont,
    label_type = "dose",
    phenotype_class = "continuous",
    alpha = 0.05,
    cutoff = 0.2,
    reliability_test = list(),
    cell_evaluation = list()
  )
  expect_equal(p_cont$family, "gaussian")
  expect_equal(as.character(p_cont$tag), paste0("dose_", seq_len(3L)))
})

test_that("ValidateScissorParams handles deprecated family and bad inputs", {
  base <- list(
    phenotype = nb_binary_phenotype(),
    label_type = "scissor",
    alpha = 0.05,
    cutoff = 0.2,
    reliability_test = list(),
    cell_evaluation = list()
  )

  expect_warning(
    p <- do.call(ValidateScissorParams, c(base, list(family = "cox"))),
    "deprecated"
  )
  expect_equal(p$phenotype_class, "survival")
  expect_equal(p$family, "cox")

  expect_error(do.call(
    ValidateScissorParams,
    c(base, list(phenotype_class = "unknown"))
  ))
  expect_error(do.call(
    ValidateScissorParams,
    c(base, list(phenotype_class = "binary", label_type = 1L))
  ))
})

test_that("ValidateScissorParams honours dots overrides", {
  p <- ValidateScissorParams(
    phenotype = nb_binary_phenotype(),
    label_type = "scissor",
    phenotype_class = "binary",
    alpha = 0.05,
    cutoff = 0.2,
    reliability_test = list(run = TRUE, n = 3L),
    cell_evaluation = list(),
    verbose = TRUE,
    seed = 7L,
    assay = "SCT"
  )

  expect_equal(p$reliability_test, list(run = TRUE, n = 3L, nfold = 10L))
  expect_true(p$verbose)
  expect_equal(p$seed, 7L)
  expect_equal(p$assay, "SCT")
})

# DoScissor - wrapper behaviour --------------------------------------------

test_that("DoScissor labels positive and negative cells from Scissor output", {
  bulk <- nb_bulk(20L, 20L, seed = 11L)
  sc <- nb_seurat(20L, 20L, seed = 12L)
  pheno <- nb_binary_phenotype(20L, seed = 13L)

  local_mocked_bindings(
    Scissor.v5.optimized = function(...) {
      list(
        Scissor_pos = colnames(sc)[1L],
        Scissor_neg = colnames(sc)[2L],
        para = list(alpha = 0.05)
      )
    },
    .package = "Scissor"
  )

  res <- DoScissor(bulk, sc, pheno, alpha = 0.05, cutoff = 0.2, verbose = FALSE, phenotype_class = "binary")

  expect_named(res, c("scRNA_data", "scissor_result", "reliability_result", "cell_evaluation"))
  expect_true(inherits(res$scRNA_data, "Seurat"))
  expect_equal(res$scRNA_data$scissor[[1L]], "Positive")
  expect_equal(res$scRNA_data$scissor[[2L]], "Negative")
  expect_equal(sum(res$scRNA_data$scissor == "Neutral"), 18L)
  expect_null(res$reliability_result)
  expect_null(res$cell_evaluation)
})

# DoScissorRelTest - early-exit branches -----------------------------------

test_that("DoScissorRelTest skips invalid parameters and single-level Y", {
  expect_warning(
    res <- DoScissorRelTest(
      scissor_res = list(Y = factor(c(0L, 1L))),
      alpha = 0.5,
      family = "binomial",
      n = 1.5,
      nfold = 10L,
      verbose = FALSE
    ),
    "Skipping reliability test"
  )
  expect_null(res)

  expect_warning(
    res <- DoScissorRelTest(
      scissor_res = list(Y = factor(rep(1L, 10L))),
      alpha = 0.5,
      family = "binomial",
      n = 2L,
      nfold = 2L,
      verbose = FALSE
    ),
    "Only one level detected"
  )
  expect_null(res)
})

# DoScissorCellEval - early-exit branches ----------------------------------

test_that("DoScissorCellEval skips missing benchmark and invalid FDR", {
  expect_warning(
    res <- DoScissorCellEval(
      benchmark_data_path = "does/not/exist.RData",
      scissor_res = list(),
      verbose = FALSE
    ),
    "does not exist"
  )
  expect_null(res)

  existing <- withr::local_tempfile(fileext = ".RData")
  file.create(existing)
  expect_warning(
    res <- DoScissorCellEval(
      benchmark_data_path = existing,
      scissor_res = list(),
      FDR_cutoff = 0,
      verbose = FALSE
    ),
    "FDR_cutoff"
  )
  expect_null(res)
})
