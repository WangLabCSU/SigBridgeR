# ---- 2. Do LP-SGL ----

#' Validate and Prepare DoLP_SGL Parameters
#'
#' @description
#' Internal helper that handles package installation checks, deprecated argument
#' warnings, input validation, and default-value resolution for [DoLP_SGL()].
#'
#' @param matched_bulk,sc_data,phenotype,label_type,phenotype_class,family
#'   Forwarded from [DoLP_SGL()].
#' @param resolution,alpha,nfold,dge_analysis Forwarded from [DoLP_SGL()].
#' @param ... Additional dots forwarded from [DoLP_SGL()].
#'
#' @return A named list with elements: `phenotype_class`, `family`, `verbose`,
#'   `seed`, `assay`, `cache_config`.
#'
#' @keywords internal
#' @family LP_SGL
ValidateLPSGLParams <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type,
  phenotype_class = c("binary", "continuous", "survival"),
  family = c("binomial", "cox", "gaussian"),
  resolution,
  alpha,
  nfold,
  dge_analysis,
  ...
) {
  # -- package checks -------------------------------------------------------
  check_installed("LPSGL", action = \(pkg, ...) {
    check_installed("pak")
    pak::pak("Exceret/LPSGL")
  })

  # -- input validation -----------------------------------------------------
  chk::chk_is(matched_bulk, c("matrix", "data.frame"))
  chk::chk_is(sc_data, c("Seurat"))
  chk::chk_list(dge_analysis)
  phenotype_class <- arg_match(phenotype_class)

  # -- handle deprecated `family` -------------------------------------------
  if (lifecycle::is_present(family)) {
    lifecycle::deprecate_warn(
      "4.0.0",
      "DoLP_SGL(family =)",
      "DoLP_SGL(phenotype_class =)"
    )
    phenotype_class <- switch(
      family,
      "binomial" = "binary",
      "cox" = "survival",
      "gaussian" = "continuous",
      Abort("Invalid family")
    )
  }
  family <- switch(
    phenotype_class,
    "binary" = "binomial",
    "survival" = "cox",
    "continuous" = "gaussian"
  )

  # -- process dots ---------------------------------------------------------
  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose")
  seed <- dots$seed %||% getFuncOption("seed")
  assay <- dots$assay %||% "RNA"

  # -- build cache config ---------------------------------------------------
  cache_config <- ScreenMethodConfig(
    method_name = "LP_SGL",
    method_version = r_pkg_version("LPSGL"),
    param = get_env_vars(exclude = c("matched_bulk", "sc_data", "phenotype")),
    phenotype_class = phenotype_class,
    label_type = label_type
  )

  get_env_vars(exclude = c("matched_bulk", "sc_data", "phenotype"))
}

#' @title Perform LP-SGL Screening Analysis
#' @description
#' Identifies phenotype-associated cell subpopulations using Lasso-Penalized
#' Sparse Group Lasso (LP-SGL) with Leiden community detection. This method
#' integrates bulk and single-cell RNA-seq data to identify cell subpopulations
#' associated with phenotypic outcomes.
#'
#' @inheritParams Screen
#' @param family `r lifecycle::badge("deprecated")` Type of regression model: "`logit`" (logistic), "`cox`" (Cox),
#' or "`linear`" (linear regression)
#' @param resolution Resolution parameter for Leiden clustering (default: `0.6`)
#' @param alpha Alpha parameter for SGL balancing L1 and L2 penalties (default: `0.5`)
#' @param nfold Number of folds for cross-validation (default: `5`)
#' @param dge_analysis List controlling differential expression analysis:
#' * `run`: Whether to run DEG analysis (default: `FALSE`)
#' * `logFC_threshold`: Log fold change threshold (default: `1`)
#' * `pval_threshold`: P-value threshold (default: `0.05`)
#' @param ... Additional arguments passed to preprocessing functions, e.g.:
#' * `verbose`: Whether to print progress messages (default: `TRUE`)
#' * `seed`: Random seed for reproducibility (default: `123L`)
#' * `assay`: Assay to use for clustering (default: `"RNA"`)
#'
#' @return A list containing:
#'
#' * **scRNA_data**: Seurat object with LP-SGL results integrated
#' * **sgl_fit**: Fitted SGL model object
#' * **cvfit**: Cross-validation results
#' * **dge_res**: Differential expression results if requested (NULL otherwise)
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
  phenotype_class = c("binary", "continuous", "survival"),
  family = lifecycle::deprecated(),
  resolution = 0.6,
  alpha = 0.5,
  nfold = 5L,
  dge_analysis = list(
    run = FALSE, # whether to run DEG analysis
    logFC_threshold = 1L,
    pval_threshold = 0.05
  ),
  ...
) {
  # -- validate & prepare all parameters -----------------------------------
  p <- exec(ValidateLPSGLParams, !!!fn_fmls())

  if (ncol(sc_data) > 5e4) {
    cli::cli_warn(
      "The number of cells in the scRNA-seq data is too large (>50k), \\
       this may result in segmentation fault"
    )
  }

  # * Start
  if (p$verbose) {
    ts_cli$cli_alert_info(cli::col_green(
      "Starting LP-SGL screen"
    ))
  }
  # * Run Leiden clustering
  leiden_results <- LPSGL::run_leiden_clustering(
    seurat_obj = sc_data,
    graph_name = paste0(p$assay, "_snn"),
    resolution = resolution,
    verbose = p$verbose,
    seed = p$seed
  )
  # * Run LP-SGL
  lpsgl_res <- LPSGL::label_cell(
    seurat_obj = sc_data,
    bulk_dataset = matched_bulk,
    phenotype = phenotype,
    cluster_membership = leiden_results,
    alpha = alpha,
    nfold = nfold,
    type = p$family,
    verbose = p$verbose,
    seed = p$seed
  )

  seurat_obj <- AddMisc(
    seurat_obj = lpsgl_res$seurat_obj,
    LP_SGL = props(p$cache_config)
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
        verbose = p$verbose
      ),
      error = function(e) {
        ts_cli$cli_alert_danger(
          cli::col_red("DEG analysis failed: ", e$message)
        )
        return(NULL)
      }
    )
  }

  if (p$verbose) {
    res_table <- table(seurat_obj$LP_SGL)
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
    scRNA_data = seurat_obj,
    sgl_fit = lpsgl_res$sgl_fit,
    cvfit = lpsgl_res$cvfit,
    dge_res = dge_res
  )
}
