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
#' @param ... Additional arguments passed to the function. Common parameters include:
#'   \describe{
#'     \item{verbose}{Logical. Whether to print verbose output (default: TRUE).}
#'   }
#'
#' @return A named list containing:
#'   \describe{
#'     \item{scRNA_data}{Modified single-cell data object with integrated screening results.}
#'   }
#'
#' @export
DoTiRank <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = "TiRank",
  phenotype_class = c("binary", "survival", "continuous"),
  tirank_params = list(
    validation_proportion = 0.15,
    sampling_thresh = 0.5,
    sampling_mode = c("smote", "downsample", "upsample", "tomeklinks"),
    top_var_genes = 2000L,
    top_gene_pairs = 1000L,
    p_value_threshold = 0.05,
    max_cutoff = 0.8,
    min_cutoff = -0.8,

    nhead = 2L,
    nhid1 = 96L,
    nhid2 = 8L,
    n_output = 32L,
    nlayers = 3L,
    n_pred = 2L,
    dropout = 0.5,
    encoder_type = c("MLP", "Transformer", "DenseNet"),
    infer_mode = c("SC", "ST"),
    n_trials = 5L,
    do_reject = TRUE,
    tolerance = 0.05
  ),
  save_path = "./TiRank_res",
  load_cache = NULL,

  ...
) {
  CheckInstalled("Exceret/rTiRank")
  # Extract additional arguments
  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose")
  seed <- dots$seed %||% getFuncOption("seed")
  assay <- dots$assay %||% "RNA"

  mode <- switch(
    phenotype_class,
    "binary" = "Classification",
    "survival" = "Cox",
    "continuous" = "Regression",
    cli::cli_abort(c(
      "x" = "`phenotype_class` must be one of 'binary', 'survival', or 'continuous'"
    ))
  )

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green("Starting TiRank Screen"))
  }

  if (is.null(load_cache)) {
    bulkExp <- rTiRank::normalize_data(matched_bulk)
    bulkClinical <- as.data.frame(phenotype)
    check_res <- rTiRank::check_bulk(
      bulkExp,
      bulkClinical,
      save_path = save_path
    )
    if (verbose) {
      ts_cli$cli_alert_info("Transfer expression profile")
    }

    st_exp_df <- rTiRank::transfer_exp_profile(sc_data)

    if (verbose) {
      ts_cli$cli_alert_info("Generate validation set")
    }

    val_res <- rTiRank::generate_val(
      check_res$bulk_exp,
      check_res$bulk_clinical,
      validation_proportion = tirank_params$validation_proportion,
      mode = mode,
      seed = seed,
    )

    if (verbose) {
      ts_cli$cli_alert_info("Perform sampling")
    }

    sampling_res <- rTiRank::perform_sampling_on_RNAseq(
      val_res$bulk_exp_train,
      val_res$bulk_clinical_train,
      mode = tirank_params$sampling_mode,
      threshold = tirank_params$sampling_thresh
    )

    if (verbose) {
      ts_cli$cli_alert_info("Compute cell-cell distance")
    }

    cell_cell_distance <- rTiRank::compute_similarity(
      seurat = seurat,
      calculate_distance = FALSE,
      parallel = FALSE,
      save_path = save_path
    )
  } else if (verbose) {
    ts_cli$cli_alert_info("Load existing data")
  }

  if (verbose) {
    ts_cli$cli_alert_info("Training model")
  }

  rTiRank::run_tirank_model(
    seurat = seurat,
    bulk_exp_train = as.data.frame(sampling_res$bulk_exp_resampled),
    bulk_clinical_train = sampling_res$bulk_clinical_resampled,
    bulk_exp_val = val_res$bulk_exp_val,
    bulk_clinical_val = val_res$bulk_clinical_val,
    save_dir = save_path,
    device = "cuda",
    gpextractor_params = list(
      top_var_genes = tirank_params$top_var_genes,
      top_gene_pairs = tirank_params$top_gene_pairs,
      p_value_threshold = tirank_params$p_value_threshold,
      max_cutoff = tirank_params$max_cutoff,
      min_cutoff = tirank_params$min_cutoff
    ),
    model_params = list(
      batch_size = tirank_params$batch_size,
      nhead = tirank_params$nhead,
      nhid1 = tirank_params$nhid1,
      nhid2 = tirank_params$nhid2,
      n_output = tirank_params$n_output,
      nlayers = tirank_params$nlayers,
      n_pred = tirank_params$n_pred,
      dropout = tirank_params$dropout,
      mode = mode,
      encoder_type = tirank_params$encoder_type,
      infer_mode = tirank_params$infer_mode,
      n_trials = tirank_params$n_trials,
      do_reject = tirank_params$do_reject,
      tolerance = tirank_params$tolerance,
      reject_mode = tirank_params$reject_mode
    ),
    sc_response_file = system.file(
      "python/Example/sc_pipeline.py",
      package = "rTiRank"
    )
  )

  # TODO: fix this read
  meta <- data.table::fread()

  modified_sc_data <- Seurat::AddMetaData(
    sc_data,
  ) %>%
    AddMisc(
      TiRank_para = tirank_params,
      TiRank_type = label_type
    )

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green("TiRank screening Done"))
  }

  # Return result in expected format
  list(
    scRNA_data = modified_sc_data,
    cell_cell_distance = cell_cell_distance
  )
}
