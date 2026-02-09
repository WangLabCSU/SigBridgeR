#' @title Perform LP-SGL Screening Analysis
#' @description
#' Identifies phenotype-associated cell subpopulations using Lasso-Penalized
#' Sparse Group Lasso (LP-SGL) with Leiden community detection. This method
#' integrates bulk and single-cell RNA-seq data to identify cell subpopulations
#' associated with phenotypic outcomes.
#'
#' @param matched_bulk Bulk expression matrix (features × samples)
#' @param sc_data Single-cell RNA-seq data (Seurat object)
#' @param phenotype Binary phenotype vector for bulk samples
#' @param label_type Character specifying phenotype label type (default: "LP_SGL")
#' @param family Type of regression model: "`logit`" (logistic), "`cox`" (Cox),
#' or "`linear`" (linear regression)
#' @param resolution Resolution parameter for Leiden clustering (default: `0.6`)
#' @param alpha Alpha parameter for SGL balancing L1 and L2 penalties (default: `0.5`)
#' @param nfold Number of folds for cross-validation (default: `5`)
#' @param dge_analysis List controlling differential expression analysis:
#' \itemize{
#' \item{run: Whether to run DEG analysis (default: `FALSE`)}
#' \item{logFC_threshold: Log fold change threshold (default: `1`)}
#' \item{pval_threshold: P-value threshold (default: `0.05`)}
#' }
#' @param ... Additional arguments passed to preprocessing functions, e.g.:
#' \itemize{
#' \item{verbose: Whether to print progress messages (default: `TRUE`)}
#' \item{seed: Random seed for reproducibility (default: `123L`)}
#' \item{assay: Assay to use for clustering (default: `"RNA"`)}
#' }
#'
#' @return A list containing:
#' \item{scRNA_data}{Seurat object with LP-SGL results integrated}
#' \item{sgl_fit}{Fitted SGL model object}
#' \item{cvfit}{Cross-validation results}
#' \item{dge_res}{Differential expression results if requested (NULL otherwise)}
#'
#' @examples
#' \dontrun{
#' # Example using simulated data
#' set.seed(123)
#'
#' # Create simulated data
#' bulk_data <- matrix(rnorm(1000*50), nrow=1000, ncol=50)
#' sc_data <- matrix(rnorm(1000*500), nrow=1000, ncol=500)
#' phenotype <- rep(c(0, 1), each=25)
#'
#' # Run LP-SGL analysis
#' results <- DoLP_SGL(
#' matched_bulk = bulk_data,
#' sc_data = sc_data,
#' phenotype = phenotype,
#' family = "logit",
#' resolution = 0.6,
#' dge_analysis = list(run = TRUE, logFC_threshold = 1, pval_threshold = 0.05)
#' )
#'
#' # Access results
#' lpsgl_seurat <- results$scRNA_data
#' sgl_model <- results$sgl_fit
#' deg_results <- results$dge_res
#' }
#'
#' @export
#' @family screen_method
#' @family LP_SGL
#' @references Li J, Zhang H, Mu B, Zuo H, Zhou K. Identifying phenotype-associated subpopulations through LP_SGL. Briefings in Bioinformatics. 2023 Nov 22;25(1):bbad424.
#'
DoLP_SGL <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = "LP_SGL",
  family = c('logit', 'cox', 'linear'),
  resolution = 0.6,
  alpha = 0.5,
  nfold = 5,
  dge_analysis = list(
    run = FALSE, # whether to run DEG analysis
    logFC_threshold = 1,
    pval_threshold = 0.05
  ),
  ...
) {
  # * Input validation
  chk::chk_is(matched_bulk, c("matrix", "data.frame"))
  chk::chk_is(sc_data, c("Seurat"))
  chk::chk_list(dge_analysis)
  family <- SigBridgeRUtils::MatchArg(
    family,
    c('logit', 'cox', 'linear'),
    NULL
  )

  # * Default params
  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose")
  seed <- dots$seed %||% getFuncOption("seed")
  assay <- dots$assay %||% "RNA"

  if (ncol(sc_data) > 5e4) {
    cli::cli_warn(
      "The number of cells in the scRNA-seq data is too large (>50k), \\
       this may result in segmentation fault"
    )
  }

  # * Start
  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green(
      "Starting LP-SGL screen"
    ))
  }
  # * Run Leiden clustering
  leiden_results <- LPSGL::run_leiden_clustering(
    seurat_obj = sc_data,
    graph_name = paste0(assay, "_snn"),
    resolution = resolution,
    verbose = verbose,
    seed = seed
  )
  # * Run LP-SGL
  lpsgl_res <- LPSGL::label_cell(
    seurat_obj = sc_data,
    bulk_dataset = matched_bulk,
    phenotype = phenotype,
    cluster_membership = leiden_results,
    alpha = alpha,
    nfold = nfold,
    type = family,
    verbose = verbose,
    seed = seed
  )

  # * Find Deferentially Expressed Genes if requested
  dge_res <- NULL
  default_dge_analysis <- list(
    run = FALSE,
    logFC_threshold = 1,
    pval_threshold = 0.05
  )
  dge_analysis <- utils::modifyList(default_dge_analysis, dge_analysis)

  if (dge_analysis$run) {
    rlang::check_installed('limma')
    dge_res <- rlang::try_fetch(
      LPSGL::perform_DEG_analysis(
        bulk_matrix = matched_bulk,
        phenotype = phenotype,
        logFC_threshold = dge_analysis$logFC_threshold,
        pval_threshold = dge_analysis$pval_threshold,
        adjust_method = 'BH',
        verbose = verbose
      ),
      error = function(e) {
        ts_cli$cli_alert_danger(
          cli::col_red("DEG analysis failed: ", e$message)
        )
        return(NULL)
      }
    )
  }

  if (verbose) {
    res_table <- table(lpsgl_res$seurat_obj$LP_SGL)
    label <- tolower(paste(label_type, names(res_table)))
    names(res_table) <- label
    msg <- purrr::imap_chr(
      res_table,
      ~ paste(.x, .y, "cells", sep = " ", collapse = ", ")
    )
    ts_cli$cli_alert_info("Identified {msg}")

    ts_cli$cli_alert_success(cli::col_green("LP-SGL screening done."))
  }

  list(
    scRNA_data = lpsgl_res$seurat_obj,
    sgl_fit = lpsgl_res$sgl_fit,
    cvfit = lpsgl_res$cvfit,
    dge_res = dge_res
  )
}
