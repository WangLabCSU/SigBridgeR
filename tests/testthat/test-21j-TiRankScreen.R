skip_if_not_installed("rTiRank")

# ---- local fixtures ---------------------------------------------------------

.tirank_args <- function(phenotype_class = "binary", ...) {
  modifyList(
    list(
      label_type = "TiRank",
      phenotype_class = phenotype_class,
      tirank_params = list(n_trials = 2L),
      save_path = "./TiRank_res",
      load_cache = NULL,
      verbose = FALSE,
      seed = 7L
    ),
    list(...)
  )
}

.tirank_params <- function() {
  list(
    validation_proportion = 0.2,
    sampling_thresh = 0.5,
    sampling_mode = "downsample",
    top_var_genes = 50L,
    top_gene_pairs = 10L,
    p_value_threshold = 0.05,
    max_cutoff = 0.8,
    min_cutoff = -0.8,
    nhead = 2L,
    nhid1 = 8L,
    nhid2 = 4L,
    n_output = 4L,
    nlayers = 1L,
    n_pred = 2L,
    dropout = 0.1,
    encoder_type = "MLP",
    infer_mode = "SC",
    n_trials = 1L,
    do_reject = TRUE,
    tolerance = 0.05,
    batch_size = 8L,
    reject_mode = "GMM"
  )
}

# `fread()` is the only file access left at the end of DoTiRank(); stub it so the
# test does not depend on the files produced by the Python backend.
.tirank_fread_mock <- function(n_cells, captured = NULL) {
  function(file) {
    if (!is.null(captured)) {
      captured$fread_file <- file
    }
    data.table::data.table(
      Reject = rep(c(TRUE, FALSE), length.out = n_cells),
      Rank_Score = seq_len(n_cells),
      Rank_Label = rep(
        c("Rank+", "Rank-", "Background"),
        length.out = n_cells
      )
    )
  }
}

# ---- ValidateTiRankParams() -------------------------------------------------

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
  args <- .tirank_args()
  args$load_cache <- NULL
  args <- c(
    args,
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

  args <- .tirank_args(load_cache = "/tmp/positional")
  expect_equal(do.call(ValidateTiRankParams, args)$load_cache, "/tmp/positional")
})

test_that("ValidateTiRankParams falls back to save_path and package options", {
  p <- do.call(ValidateTiRankParams, .tirank_args())
  expect_null(p$load_cache)
  expect_equal(p$save_cache, "./TiRank_res")

  withr::local_options(SigBridgeR.verbose = TRUE, SigBridgeR.seed = 99L)
  args <- .tirank_args()
  args$verbose <- NULL
  args$seed <- NULL
  p <- do.call(ValidateTiRankParams, args)

  expect_true(p$verbose)
  expect_equal(p$seed, 99L)
})

test_that("ValidateTiRankParams carries parameters into the cache config", {
  params <- list(n_trials = 3L, encoder_type = "Transformer")
  p <- do.call(
    ValidateTiRankParams,
    .tirank_args(label_type = "custom_label", tirank_params = params)
  )

  expect_equal(p$label_type, "custom_label")
  expect_equal(p$tirank_params, params)
  expect_equal(p$cache_config@label_type, "custom_label")
  expect_equal(p$cache_config@phenotype_class, "binary")
  expect_equal(p$cache_config@param$tirank_params, params)
  expect_equal(p$cache_config@param$save_path, "./TiRank_res")
  expect_equal(p$cache_config@method_version, utils::packageVersion("rTiRank"))
})

# ---- DoTiRank(): training path ---------------------------------------------

test_that("DoTiRank() runs the training path and adds screening metadata", {
  bulk <- nb_bulk(50L, 8L)
  pheno <- nb_binary_phenotype(8L)
  sc <- nb_seurat(50L, 6L)
  save_dir <- withr::local_tempdir()
  captured <- new.env(parent = emptyenv())

  local_mocked_bindings(
    CacheSetHere = function(...) {
      captured$cache_args <- list(...)
      save_dir
    },
    WriteCacheMeta = function(...) {
      captured$write_meta <- list(...)
      invisible(NULL)
    },
    .package = "SigBridgeR"
  )
  local_mocked_bindings(
    setup_seed = function(seed) {
      captured$seed <- seed
      invisible(NULL)
    },
    normalize_data = function(x) {
      captured$bulk <- x
      x
    },
    check_bulk = function(bulk_exp, bulk_clinical, save_path) {
      captured$check <- list(
        bulk_exp = bulk_exp,
        bulk_clinical = bulk_clinical,
        save_path = save_path
      )
      list(bulk_exp = bulk_exp, bulk_clinical = bulk_clinical)
    },
    transfer_exp_profile = function(x) {
      captured$sc <- x
      x
    },
    generate_val = function(
      bulk_exp,
      bulk_clinical,
      validation_proportion,
      mode,
      seed
    ) {
      captured$val <- list(
        validation_proportion = validation_proportion,
        mode = mode,
        seed = seed
      )
      list(
        bulk_exp_train = bulk_exp,
        bulk_clinical_train = bulk_clinical,
        bulk_exp_val = bulk_exp,
        bulk_clinical_val = bulk_clinical
      )
    },
    perform_sampling_on_RNAseq = function(
      bulk_exp_train,
      bulk_clinical_train,
      mode,
      threshold
    ) {
      captured$sampling <- list(mode = mode, threshold = threshold)
      list(
        bulk_exp_resampled = bulk_exp_train,
        bulk_clinical_resampled = bulk_clinical_train
      )
    },
    compute_similarity = function(...) {
      captured$similarity <- list(...)
      matrix(1, nrow = 2L, ncol = 2L)
    },
    run_tirank_model = function(...) {
      captured$run <- list(...)
      invisible(NULL)
    },
    .package = "rTiRank"
  )
  local_mocked_bindings(
    fread = .tirank_fread_mock(ncol(sc), captured),
    .package = "data.table"
  )

  res <- DoTiRank(
    matched_bulk = bulk,
    sc_data = sc,
    phenotype = pheno,
    phenotype_class = "binary",
    tirank_params = .tirank_params(),
    save_path = save_dir,
    verbose = FALSE,
    seed = 11L
  )

  expect_named(res, c("scRNA_data", "cell_cell_distance"))
  expect_s4_class(res$scRNA_data, "Seurat")
  expect_equal(res$cell_cell_distance, matrix(1, nrow = 2L, ncol = 2L))

  expect_equal(captured$seed, 11L)
  expect_equal(captured$cache_args$screen_method, "TiRank")
  expect_equal(captured$cache_args$mode, "save")
  expect_equal(captured$check$save_path, save_dir)
  expect_equal(captured$sc, sc)
  expect_equal(captured$val$validation_proportion, 0.2)
  expect_equal(captured$val$seed, 11L)
  expect_equal(captured$sampling$mode, "downsample")
  expect_equal(captured$sampling$threshold, 0.5)

  expect_identical(captured$run$seurat, sc)
  expect_equal(captured$run$model_params$mode, "Classification")
  expect_equal(captured$run$gpextractor_params$top_var_genes, 50L)
  expect_equal(captured$run$gpextractor_params$top_gene_pairs, 10L)
  expect_equal(captured$run$save_dir, save_dir)
  expect_null(captured$run$load_cache)

  expect_equal(
    captured$fread_file,
    file.path(save_dir, "3_Analysis/spot_predict_score.csv")
  )
  expect_equal(captured$write_meta$screen_method, "TiRank")
  expect_equal(captured$write_meta$phenotype_class, "binary")
  expect_equal(captured$write_meta$label_type, "TiRank")

  expect_true(all(
    c("TiRank_Reject", "TiRank_Rank_Score", "TiRank") %chin%
      colnames(res$scRNA_data@meta.data)
  ))
  expect_equal(
    unname(res$scRNA_data$TiRank),
    c("Positive", "Negative", "Neutral", "Positive", "Negative", "Neutral")
  )
  expect_equal(res$scRNA_data@misc$TiRank$method_name, "TiRank")
})

# ---- DoTiRank(): cache-load path -------------------------------------------

test_that("DoTiRank() loads a cache and forwards it to the model runner", {
  bulk <- nb_bulk(50L, 8L)
  pheno <- nb_binary_phenotype(8L)
  sc <- nb_seurat(50L, 6L)
  cache_dir <- withr::local_tempdir()
  captured <- new.env(parent = emptyenv())

  local_mocked_bindings(
    CacheSetHere = function(...) {
      captured$cache_args <- list(...)
      cache_dir
    },
    CheckCache = function(...) {
      captured$check <- list(...)
      invisible(TRUE)
    },
    .package = "SigBridgeR"
  )
  local_mocked_bindings(
    setup_seed = function(seed) {
      captured$seed <- seed
      invisible(NULL)
    },
    run_tirank_model = function(...) {
      captured$run <- list(...)
      # The load branch relies on a previously computed distance matrix; supply
      # one in DoTiRank()'s frame so the branch can be exercised end to end.
      assign(
        "cell_cell_distance",
        matrix(0, nrow = 1L, ncol = 1L),
        envir = parent.frame()
      )
      invisible(NULL)
    },
    .package = "rTiRank"
  )
  local_mocked_bindings(
    fread = .tirank_fread_mock(ncol(sc), captured),
    .package = "data.table"
  )

  res <- DoTiRank(
    matched_bulk = bulk,
    sc_data = sc,
    phenotype = pheno,
    phenotype_class = "binary",
    tirank_params = .tirank_params(),
    save_path = cache_dir,
    load_cache = cache_dir,
    verbose = FALSE,
    seed = 3L
  )

  expect_named(res, c("scRNA_data", "cell_cell_distance"))
  expect_equal(captured$cache_args$mode, "load")
  expect_equal(captured$check$screen_method, "TiRank")
  expect_equal(captured$run$load_cache, cache_dir)
  expect_equal(captured$run$save_dir, cache_dir)
  expect_null(captured$run$seurat)
  expect_null(captured$run$bulk_exp_train)
  expect_equal(captured$run$model_params$mode, "Classification")
  expect_equal(
    unname(res$scRNA_data$TiRank),
    c("Positive", "Negative", "Neutral", "Positive", "Negative", "Neutral")
  )
})
