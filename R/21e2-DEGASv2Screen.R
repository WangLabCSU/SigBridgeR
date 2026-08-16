# ---- 2. Do DEGASv2 ----

#' Validate and Prepare DoDEGASv2 Parameters
#'
#' @description
#' Internal helper that handles package installation checks, input validation,
#' and default-value resolution for [DoDEGASv2()].
#'
#' @param matched_bulk,phenotype,label_type,phenotype_class
#'   Forwarded from [DoDEGASv2()].
#' @param top_fraction_pos,sclab,bulk_hvg,bulk_de,sc_de,add_genes
#'   Forwarded from [DoDEGASv2()].
#' @param n_hvg,n_bulk_de,n_sc_de,padj.thresh,only.pos,min.pct,logfc.threshold
#'   Forwarded from [DoDEGASv2()].
#' @param n_st_classes,loss_type,transfer_type,model_save_dir
#'   Forwarded from [DoDEGASv2()].
#' @param lambda1,lambda2,lambda3,tot_seeds,tot_iters
#'   Forwarded from [DoDEGASv2()].
#' @param extract_embs,random_feat,random_perc,early_stopping
#'   Forwarded from [DoDEGASv2()].
#' @param ... Additional dots forwarded from [DoDEGASv2()].
#'
#' @return A named list with elements: `phenotype_class`, `loss_type`,
#'   `transfer_type`, `verbose`, `seed`.
#'
#' @keywords internal
#' @family DEGASv2
ValidateDEGASv2Params <- function(
  matched_bulk,
  phenotype,
  label_type = NULL,
  phenotype_class = c("binary", "survival", "continuous"),
  top_fraction_pos = 0.2,
  sclab = NULL,
  bulk_hvg = TRUE,
  bulk_de = TRUE,
  sc_de = TRUE,
  add_genes = NULL,
  n_hvg = 250L,
  n_bulk_de = 250L,
  n_sc_de = 200L,
  padj.thresh = 0.05,
  only.pos = FALSE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  n_st_classes = NULL,
  loss_type = c("cross_entropy", "log_neg", "rank_loss"),
  transfer_type = c("Wasserstein", "MMD"),
  model_save_dir = "DEGASv2_res",
  lambda1 = 1,
  lambda2 = 3,
  lambda3 = 3,
  tot_seeds = 10L,
  tot_iters = 300L,
  extract_embs = FALSE,
  random_feat = FALSE,
  random_perc = 0.8,
  early_stopping = FALSE,
  ...
) {
  # -- package checks -------------------------------------------------------
  check_installed("DEGASv2", action = \(pkg, ...) {
    check_installed("pak")
    pak::pak("Exceret/DEGASv2")
  })
  check_installed("dplyr")

  # -- input validation -----------------------------------------------------
  if (!is.null(label_type)) {
    chk::chk_character(label_type)
  }
  phenotype_class <- arg_match(phenotype_class)
  chk::chk_range(top_fraction_pos)
  if (!is.null(sclab)) {
    chk::chk_character(sclab)
  }
  chk::chk_matrix(matched_bulk)
  chk::chk_flag(bulk_hvg)
  chk::chk_flag(bulk_de)
  chk::chk_flag(sc_de)
  if (!is.null(add_genes)) {
    chk::chk_character(add_genes)
  }
  chk::chk_integer(
    c(n_hvg, n_bulk_de, n_sc_de),
    x_name = "n_hvg, n_bulk_de and n_sc_de"
  )
  chk::chk_range(padj.thresh)
  chk::chk_flag(only.pos)
  chk::chk_range(min.pct)
  chk::chk_number(logfc.threshold)
  if (!is.null(n_st_classes)) {
    chk::chk_integer(n_st_classes)
  }
  loss_type <- arg_match(loss_type)
  transfer_type <- arg_match(transfer_type)
  chk::chk_string(model_save_dir)
  purrr::walk(c(lambda1, lambda2, lambda3), chk::chk_number)
  chk::chk_integer(c(tot_seeds, tot_iters))
  chk::chk_flag(extract_embs)
  chk::chk_flag(random_feat)
  chk::chk_range(random_perc)
  chk::chk_flag(early_stopping)

  # -- process dots ---------------------------------------------------------
  dots <- list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose")
  seed <- dots$seed %||% getFuncOption("seed")

  # -- build cache config ---------------------------------------------------
  res <- list(
    label_type = label_type,
    phenotype_class = phenotype_class,
    top_fraction_pos = top_fraction_pos,
    sclab = sclab,
    bulk_hvg = bulk_hvg,
    bulk_de = bulk_de,
    sc_de = sc_de,
    add_genes = add_genes,
    n_hvg = n_hvg,
    n_bulk_de = n_bulk_de,
    n_sc_de = n_sc_de,
    padj.thresh = padj.thresh,
    only.pos = only.pos,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    n_st_classes = n_st_classes,
    loss_type = loss_type,
    transfer_type = transfer_type,
    model_save_dir = file.path(getwd(), model_save_dir),
    lambda1 = lambda1,
    lambda2 = lambda2,
    lambda3 = lambda3,
    tot_seeds = tot_seeds,
    tot_iters = tot_iters,
    extract_embs = extract_embs,
    random_feat = random_feat,
    random_perc = random_perc,
    early_stopping = early_stopping,
    verbose = verbose,
    seed = seed,
    dots = dots
  )
  cache_config <- ScreenMethodConfig(
    method_name = "DEGASv2",
    param = res,
    phenotype_class = phenotype_class,
    label_type = label_type
  )
  # -- return validated parameters ------------------------------------------
  c(res, list(cache_config = cache_config))
}

#' @title Screen Function Template
#'
#' @param matched_bulk Matrix or data frame of preprocessed bulk RNA-seq expression
#'        data (genes x samples). Column names must match names/IDs in `phenotype`.
#' @param sc_data A matrix/Matrix (genes x cells) or a Seurat object containing scRNA-seq data to be screened.
#' @param phenotype Phenotype data, either:
#'        - Named vector (names match `matched_bulk` columns), or
#'        - Patient survival Data frame with row names matching `matched_bulk` columns, colnames named "time" and "status"
#' @param label_type Character specifying phenotype label type
#' @param phenotype_class Type of phenotypic outcome (must be consistent with input data):
#'        - `"binary"`: Binary traits (e.g., case/control)
#'        - `"continuous"`: Continuous measurements
#'        - `"survival"`: Survival infomation
#' @param top_fraction_pos Proportion of cells to be labelled as `"Positive"`
#' @inheritParams DEGASv2::DEGAS_preprocessing
#' @inheritParams DEGASv2::run_DEGAS_SCST
#' @inheritParams Seurat::FindAllMarkers
#'
#' @param ... Additional arguments passed to the function. Common parameters include:
#'    - verbose: Logical. Whether to print verbose output (default: TRUE).
#'    - seed: Integer. Random seed for reproducibility (default: 123).
#'    - Other parameters are passed to [Seurat::FindAllMarkers()]
#'
#'
#' @return A named list containing:
#'   \describe{
#'     \item{scRNA_data}{Modified single-cell data object with integrated screening results.}
#'   }
#'
#' @export
DoDEGASv2 <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = NULL,
  phenotype_class = c("binary", "survival", "continuous"),
  top_fraction_pos = 0.2,
  sclab = NULL,
  bulk_hvg = TRUE,
  bulk_de = TRUE,
  sc_de = TRUE,
  add_genes = NULL,
  n_hvg = 250L,
  n_bulk_de = 250L,
  n_sc_de = 200L,
  padj.thresh = 0.05,
  only.pos = FALSE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  n_st_classes = length(unique(sc_data$seurat_clusters)),
  loss_type = c("cross_entropy", "log_neg", "rank_loss"),
  transfer_type = c("Wasserstein", "MMD"),
  model_save_dir = "DEGASv2_res",
  lambda1 = 1,
  lambda2 = 3,
  lambda3 = 3,
  tot_seeds = 10L,
  tot_iters = 300L,
  extract_embs = FALSE,
  random_feat = FALSE,
  random_perc = 0.8,
  early_stopping = FALSE,
  ...
) {
  # -- validate & prepare all parameters -----------------------------------
  p <- ValidateDEGASv2Params(
    matched_bulk = matched_bulk,
    phenotype = phenotype,
    label_type = label_type,
    phenotype_class = phenotype_class,
    top_fraction_pos = top_fraction_pos,
    sclab = sclab,
    bulk_hvg = bulk_hvg,
    bulk_de = bulk_de,
    sc_de = sc_de,
    add_genes = add_genes,
    n_hvg = n_hvg,
    n_bulk_de = n_bulk_de,
    n_sc_de = n_sc_de,
    padj.thresh = padj.thresh,
    only.pos = only.pos,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    n_st_classes = n_st_classes,
    loss_type = loss_type,
    transfer_type = transfer_type,
    model_save_dir = model_save_dir,
    lambda1 = lambda1,
    lambda2 = lambda2,
    lambda3 = lambda3,
    tot_seeds = tot_seeds,
    tot_iters = tot_iters,
    extract_embs = extract_embs,
    random_feat = random_feat,
    random_perc = random_perc,
    early_stopping = early_stopping,
    ...
  )

  set.seed(p$seed)

  if (p$verbose) {
    ts_cli$cli_alert_info(cli::col_green("Start DEGASv2 screening"))
  }

  umap_df <- umap_coordinate(obj = sc_data)

  preprocessing_model_type <- if (p$phenotype_class == "survival") {
    "survival"
  } else {
    "category"
  }

  data <- DEGAS_preprocessing(
    scst_list = sc_data,
    sclab = sclab,
    patdata = matched_bulk,
    phenotype = phenotype,
    bulk_hvg = bulk_hvg,
    bulk_de = bulk_de,
    sc_de = sc_de,
    add_genes = add_genes,
    n_hvg = n_hvg,
    n_bulk_de = n_bulk_de,
    n_sc_de = n_sc_de,
    padj.thresh = padj.thresh,
    model_type = preprocessing_model_type,
    ...
  )

  model_type <- get_model_type(phenotype = phenotype, sclab = sclab)

  if (p$verbose) {
    ts_cli$cli_alert_info("Training DEGASv2 model")
  }

  # purrr::iwalk(data, ~ cli::cli_alert_info("dim of {.y}: {dim(.x)}"))

  degas_sc_results <- DEGASv2::run_DEGAS_SCST(
    data_list = data,
    model_type = model_type,
    data_name = label_type,
    loss_type = p$loss_type,
    transfer_type = p$transfer_type,
    model_save_dir = model_save_dir,
    lambda1 = lambda1,
    lambda2 = lambda2,
    lambda3 = lambda3,
    tot_seeds = tot_seeds,
    tot_iters = tot_iters,
    extract_embs = extract_embs,
    random_feat = random_feat,
    random_perc = random_perc,
    early_stopping = early_stopping,
    n_st_classes = n_st_classes
  )
  hazard_df <- as.data.frame(degas_sc_results)
  colnames(hazard_df) <- paste0("DEGASv2_", colnames(hazard_df))

  top_fraction <- stats::quantile(
    hazard_df$DEGASv2_hazard,
    1 - top_fraction_pos
  )
  hazard_df <- dplyr::mutate(
    hazard_df,
    DEGASv2_index = NULL,
    DEGASv2_pos_thresh = top_fraction,
    DEGASv2 = ifelse(
      hazard_df$DEGASv2_hazard > top_fraction,
      "Positive",
      "Other"
    )
  )

  sc_data <- Seurat::AddMetaData(
    object = sc_data,
    metadata = hazard_df
  )
  sc_data <- AddMisc(
    seurat_obj = sc_data,
    DEGASv2 = props(p$cache_config),
    cover = FALSE
  )

  if (p$verbose) {
    ts_cli$cli_alert_info(cli::col_green("DEGASv2 screening done"))
  }

  list(
    scRNA_data = sc_data
  )
}

DEGAS_preprocessing <- function(
  scst_list, # a seurat obj
  patdata, # bulk RNA-seq data
  phenotype,
  sclab = NULL, # seurat meta
  bulk_hvg = TRUE,
  bulk_de = TRUE,
  sc_de = TRUE,
  add_genes = NULL,
  n_hvg = 250L,
  n_bulk_de = 250L,
  n_sc_de = 200L,
  padj.thresh = 0.05,
  model_type = c("category", "survival"),
  only.pos = FALSE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = TRUE,
  ...
) {
  gene_list <- select_genes2(
    scdata = scst_list,
    # sclab = sclab,
    patdata = patdata,
    phenotype = phenotype,
    add_genes = add_genes,
    bulk_hvg = bulk_hvg,
    bulk_de = bulk_de,
    sc_de = sc_de,
    n_hvg = n_hvg,
    n_bulk_de = n_bulk_de,
    n_sc_de = n_sc_de,
    padj.thresh = padj.thresh,
    model_type = model_type,
    verbose = verbose,
    ...
  )

  if (verbose) {
    ts_cli$cli_alert_info(
      "Normalizing data"
    )
  }

  norm_out <- normalize_counts_with_selected_genes(
    bulk_dataset = patdata,
    scst_list = scst_list,
    gene_list = gene_list # a chr vec
  )

  # clean phenotype
  if (model_type != "survival") {
    phenotype <- as.factor(phenotype)
    phenotype <- as.integer(phenotype) - 1
  }

  if (!is.null(sclab)) {
    sclab <- as.integer(as.factor(sclab)) - 1
    message("sclab: ", toString(unique(sclab)), "\n")
  }

  list(
    patDat = norm_out$patDat, # transposed
    phenotype = phenotype,
    scstDat = norm_out$scstDat, # transposed
    scstName = norm_out$scstName,
    sclab = sclab
  )
}


select_genes2 <- function(
  scdata,
  # sclab,
  patdata,
  phenotype,
  add_genes = NULL,
  bulk_hvg = TRUE,
  bulk_de = TRUE,
  sc_de = TRUE,
  n_hvg = 250L,
  n_bulk_de = 250L,
  n_sc_de = 200L,
  padj_thresh = 0.05,
  model_type = c("category", "survival"),
  verbose = TRUE,
  ...
) {
  genes <- c()

  # * Bulk HVG
  if (isTRUE(bulk_hvg)) {
    if (verbose) {
      ts_cli$cli_alert_info("Finding highly variable genes in bulk data")
    }

    # cpp fn
    bulk_hvgs <- get_bulk_hvg(
      patdata = patdata,
      common_genes = rownames(patdata),
      n_hvg = n_hvg
    )
    genes <- union(genes, bulk_hvgs)
  }

  # * Bulk DE
  if (isTRUE(bulk_de)) {
    if (!IsCountsMatrix(patdata, verbose = FALSE)) {
      cli::cli_alert_warning(cli::col_yellow(
        "Bulk data is not a counts matrix, skipping bulk differential expression analysis"
      ))
    } else {
      if (verbose) {
        ts_cli$cli_alert_info(
          "Finding differentially expressed genes in bulk data"
        )
      }
      bulk_de_genes <- get_bulk_deg(
        patdata = patdata,
        phenotype = phenotype,
        model_type = model_type
      )
      genes <- union(genes, bulk_de_genes)
    }
  }

  # * SC DE
  if (isTRUE(sc_de)) {
    if (verbose) {
      ts_cli$cli_alert_info(
        "Finding differentially expressed genes in single-cell data"
      )
    }
    sc_de_genes <- get_sc_deg(
      obj = scdata,
      only.pos = FALSE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      n_sc_de = n_sc_de,
      verbose = verbose,
      ...
    )
    genes <- union(genes, sc_de_genes)
  }

  # * whether to add individual gene list
  if (!is.null(add_genes)) {
    add_genes <- as.character(add_genes)
    genes <- union(genes, intersect(add_genes, rownames(patdata)))
  }
  genes <- unique(genes[!is.na(genes)])

  genes
}

get_bulk_deg <- function(
  patdata, # bulk counts matrix
  phenotype,
  model_type = c("category", "survival"),
  padj_thresh = 0.05,
  n_bulk_de = 250L
) {
  pat_lab <- if (model_type == "survival") {
    surv_mid <- stats::median(phenotype$time)
    (phenotype$time < surv_mid) * phenotype$status + 1L
  } else {
    phenotype + 1L
  }
  pat_lab <- as.factor(pat_lab)

  high_diff_genes <- tryCatch(
    {
      patdata[patdata < 0L] <- 0L

      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = patdata,
        colData = data.frame(id = colnames(patdata), label = pat_lab),
        design = ~label
      )
      dds <- DESeq2::DESeq(dds)
      res <- stats::na.omit(DESeq2::results(dds))
      res <- res[order(res$padj), ]
      res <- res[res$padj < padj_thresh & !is.na(res$padj), ]

      rownames(res)[seq_len(min(n_bulk_de, nrow(res)))]
    },
    error = function(e) {
      warning("Bulk DEG selection failed and will be skipped: ", e$message)
      return(NULL)
    }
  )
  high_diff_genes
}

get_sc_deg <- function(
  obj,
  only.pos = FALSE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  n_sc_de = 200L,
  ...
) {
  check_installed(c("dplyr", "DESeq2"))
  p_val_adj <- avg_log2FC <- gene <- NULL # suppress checking NOTE

  sc_markers <- Seurat::FindAllMarkers(
    obj,
    only.pos = only.pos,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    ...
  )
  stats::na.omit(sc_markers) |>
    dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC)) |>
    dplyr::top_n(n_sc_de) |>
    dplyr::pull(gene)
}

normalize_counts_with_selected_genes <- function(
  bulk_dataset,
  scst_list, # a single cell expr counts matrix
  gene_list # a chr vec
) {
  # normalize bulk
  bulk_dataset <- bulk_dataset[rownames(bulk_dataset) %chin% gene_list, ]
  st_counts <- SeuratObject::LayerData(
    scst_list,
    assay = "RNA",
    layer = "counts"
  )
  st_counts <- st_counts[rownames(st_counts) %chin% gene_list, ] # single cell expr counts

  common_genes <- intersect(rownames(bulk_dataset), rownames(scst_list))
  bulk_dataset <- bulk_dataset[common_genes, ]
  st_counts <- st_counts[common_genes, ]
  st_counts <- as.matrix(st_counts) # S4Matrix -> matrix

  pat_dat <- preprocessCounts_ptr(bulk_dataset) # cpp fn
  colnames(pat_dat) <- rownames(bulk_dataset)
  rownames(pat_dat) <- colnames(bulk_dataset)

  scst_expr_mat <- preprocessCounts_ptr(st_counts)
  colnames(scst_expr_mat) <- rownames(st_counts)
  rownames(scst_expr_mat) <- colnames(st_counts)

  scst_names <- rep("Dataset1", ncol(st_counts)) # number of genes

  return(list(
    patDat = pat_dat, # transposed
    scstDat = scst_expr_mat, # transposed
    scstName = scst_names
  ))
}


umap_coordinate <- function(obj) {
  umap_df <- Seurat::Embeddings(obj, "umap") |> as.data.frame()
  colnames(umap_df) <- c("UMAP_1", "UMAP_2")
  umap_df$cell <- rownames(umap_df)

  return(umap_df)
}

get_model_type <- function(phenotype, sclab = NULL) {
  bulk_level <- if (is_2d(phenotype)) {
    "Cox"
  } else {
    "Class"
  }

  sc_level <- if (!is.null(sclab)) {
    "Class"
  } else {
    "Blank"
  }

  paste0(sc_level, bulk_level)
}
