#' @title Perform TiRank Screening Analysis
#'
#' @description
#' Integrates single-cell and bulk RNA-seq data using the TiRank deep learning
#' framework to identify phenotype-associated cells. TiRank employs a
#' rank-based approach with neural network models to score and classify cells
#' based on their association with the phenotype of interest.
#'
#' @param matched_bulk Matrix or data frame of preprocessed bulk RNA-seq expression
#'        data (genes x samples). Column names must match names/IDs in `phenotype`.
#' @param sc_data A matrix/Matrix (genes x cells) or a Seurat object containing
#'        scRNA-seq data to be screened.
#' @param phenotype Phenotype data, either:
#'        - Named vector (names match `matched_bulk` columns), or
#'        - Patient survival data frame with row names matching `matched_bulk`
#'          columns, colnames named "time" and "status"
#' @param label_type Character specifying phenotype label type (default: `"TiRank"`).
#' @param phenotype_class Type of phenotypic outcome (must be consistent with input data):
#'        - `"binary"`: Binary traits (e.g., case/control)
#'        - `"continuous"`: Continuous measurements
#'        - `"survival"`: Survival information
#' @param tirank_params List of TiRank algorithm parameters:
#'   \describe{
#'     \strong{Data preprocessing:}
#'     \item{validation_proportion}{Proportion of bulk data held out for validation (default: `0.15`).}
#'     \item{sampling_thresh}{Threshold for resampling strategy (default: `0.5`).}
#'     \item{sampling_mode}{Resampling method: `"smote"`, `"downsample"`, `"upsample"`, or `"tomeklinks"`.}
#'     \item{top_var_genes}{Number of top variable genes to select (default: `2000L`).}
#'     \item{top_gene_pairs}{Number of top gene pairs for ranking (default: `1000L`).}
#'     \item{p_value_threshold}{P-value threshold for feature selection (default: `0.05`).}
#'     \item{max_cutoff}{Upper cutoff for correlation filtering (default: `0.8`).}
#'     \item{min_cutoff}{Lower cutoff for correlation filtering (default: `-0.8`).}
#'     \strong{Neural network architecture:}
#'     \item{nhead}{Number of attention heads (default: `2L`).}
#'     \item{nhid1}{Size of first hidden layer (default: `96L`).}
#'     \item{nhid2}{Size of second hidden layer (default: `8L`).}
#'     \item{n_output}{Output dimension of encoder (default: `32L`).}
#'     \item{nlayers}{Number of encoder layers (default: `3L`).}
#'     \item{n_pred}{Number of prediction heads (default: `2L`).}
#'     \item{dropout}{Dropout rate for regularization (default: `0.5`).}
#'     \item{encoder_type}{Encoder architecture: `"MLP"`, `"Transformer"`, or `"DenseNet"`.}
#'     \item{infer_mode}{Inference mode: `"SC"` (single-cell) or `"ST"` (spatial transcriptomics).}
#'     \strong{Training:}
#'     \item{n_trials}{Number of repeated training trials (default: `5L`).}
#'     \item{do_reject}{Whether to apply rejection criteria to uncertain predictions (default: `TRUE`).}
#'     \item{tolerance}{Tolerance threshold for rejection (default: `0.05`).}
#'   }
#' @param save_path Directory path for saving intermediate and final results
#'        (default: `"./TiRank_res"`).
#' @param load_cache Optional path to cached data for re-running analysis
#'        without recomputing expensive steps (default: `NULL`).
#' @param ... Additional arguments passed to the function. Common parameters include:
#'   \describe{
#'     \item{verbose}{Logical. Whether to print verbose output (default: TRUE).}
#'     \item{seed}{Integer. Random seed for reproducibility.}
#'     \item{assay}{Character. Name of assay to use from Seurat object (default: `"RNA"`).}
#'   }
#'
#' @return A named list containing:
#'   \describe{
#'     \item{scRNA_data}{Modified single-cell data object with integrated screening
#'       results added as metadata, including \code{TiRank_Reject}, \code{TiRank_Rank_Score},
#'       and \code{TiRank} (Positive/Neutral/Negative) columns, plus \code{TiRank_para}
#'       and \code{TiRank_type} stored in misc slot.}
#'     \item{cell_cell_distance}{Computed cell-cell similarity/distance matrix.}
#'   }
#'
#' @details
#' The TiRank screening workflow consists of the following steps:
#' \enumerate{
#'   \item \strong{Data preprocessing:} Normalizes bulk expression data and validates
#'         compatibility with phenotype information.
#'   \item \strong{Expression transfer:} Transfers the single-cell expression profile
#'         into the TiRank-compatible format.
#'   \item \strong{Validation set generation:} Splits bulk data into training and
#'         validation sets according to \code{validation_proportion}.
#'   \item \strong{Resampling:} Applies the specified sampling strategy to address
#'         class imbalance in bulk data.
#'   \item \strong{Cell-cell similarity:} Computes cell-cell distance or similarity
#'         matrix from single-cell data.
#'   \item \strong{Model training:} Trains a TiRank neural network (MLP, Transformer,
#'         or DenseNet encoder) using bulk expression and phenotype.
#'   \item \strong{Screening:} Applies the trained model to score each cell, producing
#'         rank-based predictions with rejection handling.
#'   \item \strong{Label assignment:} Classifies cells as \code{"Positive"},
#'         \code{"Neutral"}, or \code{"Negative"} based on rank scores.
#' }
#'
#' @family TiRank
#' @family screen_method
#' @export
#'
#' @examples
#' \dontrun{
#' # Binary classification example
#' result <- DoTiRank(
#'   matched_bulk = bulk_matrix,
#'   sc_data = seurat_obj,
#'   phenotype = bulk_phenotype,
#'   phenotype_class = "binary",
#'   save_path = "./TiRank_res"
#' )
#'
#' # Survival analysis example
#' result <- DoTiRank(
#'   matched_bulk = bulk_matrix,
#'   sc_data = seurat_obj,
#'   phenotype = survival_data,
#'   phenotype_class = "survival",
#'   tirank_params = list(
#'     encoder_type = "Transformer",
#'     n_trials = 3L
#'   ),
#'   save_path = "./TiRank_survival_res"
#' )
#'
#' # Access results
#' modified_seurat <- result$scRNA_data
#' head(modified_seurat[[]])
#' }
#'
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
    chk::chk_not_null(save_path)
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
      seurat = sc_data,
      calculate_distance = FALSE,
      parallel = FALSE,
      save_path = save_path
    )
    rm(cell_cell_distance)
  } else if (verbose) {
    ts_cli$cli_alert_info("Load existing data")
  }

  if (verbose) {
    ts_cli$cli_alert_info("Training model")
  }

  rTiRank::run_tirank_model(
    seurat = sc_data,
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
  meta <- data.table::fread(file.path(
    save_path,
    "3_Analysis/spot_predict_score.csv"
  ))
  meta_to_add <- meta[,
    TiRank := data.table::fcase(
      Rank_Label == "Background" , "Neutral"  ,
      Rank_Label == "Rank-"      , "Negative" ,
      Rank_Label == "Rank+"      , "Positive"
    )
  ][,
    .(Reject, Rank_Score, TiRank)
  ]

  data.table::setnames(
    meta_to_add,
    old = names(meta_to_add)[1:2],
    new = paste0("TiRank_", names(meta_to_add)[1:2])
  )

  modified_sc_data <- Seurat::AddMetaData(
    sc_data,
    meta_to_add
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
