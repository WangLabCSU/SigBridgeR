skip_if_not_installed("DEGASv2")
skip_if_not_installed("dplyr")
skip_if_not_installed("DESeq2")

.degasv2_args <- function(phenotype_class = "binary") {
  list(
    matched_bulk = nb_bulk(),
    phenotype = nb_binary_phenotype(),
    label_type = "DEGASv2",
    phenotype_class = phenotype_class,
    top_fraction_pos = 0.2,
    sclab = NULL,
    bulk_hvg = TRUE,
    bulk_de = FALSE,
    sc_de = FALSE,
    add_genes = NULL,
    n_hvg = 50L,
    n_bulk_de = 50L,
    n_sc_de = 50L,
    padj.thresh = 0.05,
    only.pos = FALSE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    n_st_classes = 2L,
    loss_type = "cross_entropy",
    transfer_type = "Wasserstein",
    model_save_dir = "DEGASv2_res",
    lambda1 = 1L,
    lambda2 = 3L,
    lambda3 = 3L,
    tot_seeds = 2L,
    tot_iters = 5L,
    extract_embs = FALSE,
    random_feat = FALSE,
    random_perc = 0.8,
    early_stopping = FALSE,
    verbose = FALSE,
    seed = 7L
  )
}

test_that("ValidateDEGASv2Params resolves defaults", {
  p <- do.call(ValidateDEGASv2Params, .degasv2_args())

  expect_equal(p$phenotype_class, "binary")
  expect_equal(p$loss_type, "cross_entropy")
  expect_equal(p$transfer_type, "Wasserstein")
  expect_equal(p$model_save_dir, file.path(getwd(), "DEGASv2_res"))
  expect_equal(p$seed, 7L)
  expect_false(p$verbose)
  expect_true(S7::S7_inherits(p$cache_config, ScreenMethodConfig))
  expect_equal(p$cache_config@method_name, "DEGASv2")
})

test_that("ValidateDEGASv2Params validates algorithm choices", {
  expect_equal(
    do.call(ValidateDEGASv2Params, c(.degasv2_args(), list(phenotype_class = "survival")))$phenotype_class,
    "survival"
  )
  expect_error(do.call(
    ValidateDEGASv2Params,
    within(.degasv2_args(), loss_type <- "invalid")
  ))
  expect_error(do.call(
    ValidateDEGASv2Params,
    within(.degasv2_args(), transfer_type <- "invalid")
  ))
  expect_error(do.call(
    ValidateDEGASv2Params,
    within(.degasv2_args(), matched_bulk <- as.data.frame(nb_bulk(10L, 10L)))
  ))
  expect_error(do.call(
    ValidateDEGASv2Params,
    within(.degasv2_args(), top_fraction_pos <- 2)
  ))
})

test_that("get_model_type derives model type from phenotype shape", {
  expect_equal(get_model_type(nb_binary_phenotype()), "BlankClass")
  expect_equal(get_model_type(nb_survival_phenotype()), "BlankCox")
  expect_equal(get_model_type(nb_binary_phenotype(), sclab = c("A", "B")), "ClassClass")
  expect_equal(get_model_type(nb_survival_phenotype(), sclab = "A"), "ClassCox")
})

test_that("umap_coordinate returns a data frame with cell names", {
  sc <- nb_seurat(20L, 20L)
  embeddings <- matrix(
    stats::rnorm(40L),
    nrow = 20L,
    dimnames = list(colnames(sc), c("UMAP_1", "UMAP_2"))
  )
  sc[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings,
    key = "UMAP_"
  )

  df <- umap_coordinate(sc)

  expect_s3_class(df, "data.frame")
  expect_equal(names(df), c("UMAP_1", "UMAP_2", "cell"))
  expect_equal(df$cell, colnames(sc))
})

test_that("normalize_counts_with_selected_genes transposes both datasets", {
  bulk <- nb_bulk(50L, 60L)
  sc <- nb_seurat(50L, 70L)
  genes <- rownames(bulk)[1L:30L]

  out <- normalize_counts_with_selected_genes(bulk, sc, genes)

  expect_equal(dim(out$patDat), c(60L, 30L))
  expect_equal(dim(out$scstDat), c(70L, 30L))
  expect_equal(out$scstName, rep("Dataset1", 70L))
})

test_that("select_genes2 unions bulk HVGs and user-supplied genes", {
  bulk <- nb_bulk(50L, 50L)
  sc <- nb_seurat(50L, 50L)

  genes <- select_genes2(
    scdata = sc,
    patdata = bulk,
    phenotype = nb_binary_phenotype(50L),
    add_genes = c(rownames(bulk)[1L], "NOT_A_GENE"),
    bulk_hvg = TRUE,
    bulk_de = FALSE,
    sc_de = FALSE,
    n_hvg = 20L,
    verbose = FALSE
  )

  expect_true(is.character(genes))
  expect_true(length(genes) > 0L)
  expect_true("G1" %chin% genes)
  expect_false("NOT_A_GENE" %chin% genes)
  expect_setequal(genes, unique(genes))
})
