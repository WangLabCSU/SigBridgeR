# ---- 2. Do scPAS ----

#' Validate and Prepare DoscPAS Parameters
#'
#' @description
#' Internal helper that handles package installation checks, deprecated argument
#' warnings, input validation, and default-value resolution for [DoscPAS()].
#'
#' @param phenotype,label_type,phenotype_class,assay
#'   Forwarded from [DoscPAS()].
#' @param imputation,imputation_method,nfeature,alpha,cutoff,network_class
#'   Forwarded from [DoscPAS()].
#' @param family Deprecated argument forwarded from [DoscPAS()].
#' @param permutation_times,FDR_threshold,independent Forwarded from [DoscPAS()].
#' @param ... Additional dots forwarded from [DoscPAS()].
#'
#' @return A named list with elements: `label_type`, `family`, `imputation_method`,
#'   `network_class`, `verbose`, `seed`.
#'
#' @keywords internal
#' @family scPAS
ValidatescPASParams <- function(
  phenotype,
  label_type,
  phenotype_class,
  assay,
  imputation,
  imputation_method = c("KNN", "ALRA"),
  nfeature,
  alpha,
  cutoff,
  network_class = c("SC", "bulk"),
  family = lifecycle::deprecated(),
  permutation_times,
  FDR_threshold,
  independent,
  ...
) {
  # -- package checks -------------------------------------------------------
  check_installed("dplyr")
  check_installed("scPAS", action = \(pkg, ...) {
    check_installed("pak")
    pak::pak("Exceret/scPAS")
  })

  # -- input validation -----------------------------------------------------
  chk::chk_character(label_type)
  chk::chk_flag(imputation)
  chk::chk_flag(independent)
  chk::chk_number(permutation_times)
  chk::chk_number(FDR_threshold)
  chk::chk_number(nfeature)
  chk::chk_range(cutoff)
  chk::chk_length(label_type)

  if (imputation) {
    imputation_method <- arg_match(imputation_method)
  }
  if (!is.null(alpha)) {
    chk::chk_range(alpha)
  }
  network_class <- arg_match(network_class)

  # -- handle deprecated `family` -------------------------------------------
  if (lifecycle::is_present(family)) {
    lifecycle::deprecate_warn(
      "4.0.0",
      "DoscPAS(family = )",
      "DoscPAS(phenotype_class = )"
    )
    phenotype_class <- switch(
      family,
      "cox" = "survival",
      "gaussian" = "continuous",
      "binomial" = "binary",
      Abort("Invalid family: {.val {family}}")
    )
  }
  family <- switch(
    phenotype_class,
    "binary" = "binomial",
    "continuous" = "gaussian",
    "survival" = "cox",
    Abort("Invalid phenotype_class: {.val {phenotype_class}}")
  )
  label_type <- paste0(label_type, c("_Positive", "_Negative"))

  # -- process dots ---------------------------------------------------------
  dots <- list2(...)
  verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption("verbose")
  seed <- dots$seed %||% SigBridgeRUtils::getFuncOption("seed")
  load_cache <- dots$load_cache
  save_cache <- dots$save_cache

  res <- list(
    label_type = label_type,
    phenotype_class = phenotype_class,
    family = family,
    assay = assay,
    imputation = imputation,
    imputation_method = imputation_method,
    nfeature = nfeature,
    alpha = alpha,
    cutoff = cutoff,
    network_class = network_class,
    permutation_times = permutation_times,
    FDR_threshold = FDR_threshold,
    independent = independent,
    verbose = verbose,
    seed = seed,
    load_cache = load_cache,
    save_cache = save_cache,
    dots = dots
  )
  cache_config <- ScreenMethodConfig(
    method_name = "scPAS",
    param = res,
    phenotype_class = phenotype_class,
    label_type = label_type
  )

  c(res, list(cache_config = cache_config))
}

#' @title Perform scPAS Screening Analysis
#' @description
#' This function performs scPAS screening analysis by integrating bulk and single-cell RNA-seq data.
#' It includes data filtering steps and wraps the core scPAS::scPAS function.
#'
#' @inheritParams Screen
#' @param assay Assay to use from sc_data (default: 'RNA')
#' @param imputation Logical, whether to perform imputation (default: FALSE)
#' @param imputation_method Character. Name of alternative method for imputation. (options: "KNN", "ALRA")
#' @param nfeature Number of features to select (default: 3000, indicating that the top 3000 highly variable genes are selected for model training
#' @param alpha Numeric. Significance threshold. Parameter used to balance the effect of the l1 norm and the network-based penalties. It can be a number or a searching vector. If alpha = NULL, a default searching vector is used. The range of alpha is in `[0,1]`. A larger alpha lays more emphasis on the l1 norm. (default: 0.01)
#' @param cutoff Numeric. Cutoff value for selecting the optimal alpha value when alpha = NULL. (default: 0.2)
#' @param network_class Network class to use (default: 'SC', indicating gene-gene similarity networks derived from single-cell data. The other one is 'bulk'.)
#' @param family `r lifecycle::badge("deprecated")` Model family for analysis (options: "cox", "gaussian", "binomial")
#' @param permutation_times Number of permutations to perform (default: 2000)
#' @param FDR_threshold Numeric. FDR value threshold for identifying phenotype-associated cells (default: 0.05)
#' @param independent Logical. The background distribution of risk scores is constructed independently of each cell. (default: TRUE)
#' @param ... Additional arguments. Currently supports:
#'    - `verbose`: Logical indicating whether to print progress messages. Defaults to `TRUE`.
#'    - `seed`: For reproducibility, default is `123L`
#'
#' @return A Seurat object from scPAS analysis
#'
#' @references
#' Xie A, Wang H, Zhao J, Wang Z, Xu J, Xu Y. scPAS: single-cell phenotype-associated subpopulation identifier. Briefings in Bioinformatics. 2024 Nov 22;26(1):bbae655.
#'
#' @section LICENSE:
#' Licensed under the GNU General Public License version 3 (GPL-3.0).
#' A copy of the license is available at <https://www.gnu.org/licenses/gpl-3.0.en.html>.
#'
#' @export
#'
#' @family screen_method
#' @family scPAS
#'
DoscPAS <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = "scPAS",
  phenotype_class = c("binary", "continuous", "survival"),
  assay = 'RNA',
  imputation = FALSE,
  imputation_method = c("KNN", "ALRA"),
  nfeature = 3000L,
  alpha = c(0.01, NULL),
  cutoff = 0.2,
  network_class = c("SC", "bulk"),
  family = lifecycle::deprecated(),
  permutation_times = 2000L,
  FDR_threshold = 0.05,
  independent = TRUE,
  ...
) {
  # -- validate & prepare all parameters -----------------------------------
  p <- ValidatescPASParams(
    phenotype = phenotype,
    label_type = label_type,
    phenotype_class = phenotype_class,
    assay = assay,
    imputation = imputation,
    imputation_method = imputation_method,
    nfeature = nfeature,
    alpha = alpha,
    cutoff = cutoff,
    network_class = network_class,
    family = family,
    permutation_times = permutation_times,
    FDR_threshold = FDR_threshold,
    independent = independent,
    ...
  )

  set.seed(p$seed)

  if (p$verbose) {
    ts_cli$cli_alert_info(cli::col_green("Start scPAS screening."))
  }

  # Set default assay
  Seurat::DefaultAssay(sc_data) <- assay

  # * Step 0: Common gene identification with optimized filtering
  bulk_dataset <- as.matrix(matched_bulk)
  common_genes <- scPAS::identify_common_genes(
    bulk_dataset = bulk_dataset,
    sc_dataset = sc_data,
    nfeature = nfeature,
    verbose = p$verbose
  )

  # Step 1: Quantile normalization with matrix optimization
  if (p$verbose) {
    ts_cli$cli_alert_info("Quantile normalizing bulk data")
  }

  Expression_bulk <- SigBridgeRUtils::normalize.quantiles(
    bulk_dataset[common_genes, ],
    keep.names = TRUE
  )

  # Step 2: Single-cell expression processing
  if (imputation) {
    sc_data <- scPAS::imputation2(
      obj = sc_data,
      assay = assay,
      method = p$imputation_method,
      verbose = p$verbose
    )
    assay <- Seurat::DefaultAssay(sc_data)
  }
  if (p$verbose) {
    ts_cli$cli_alert_info(
      "Extracting single-cell expression profiles"
    )
  }

  sc_exprs <- SeuratObject::LayerData(
    sc_data,
    assay = assay,
    layer = 'data'
  ) # Get expression data from Seurat
  Expression_cell <- sc_exprs[common_genes, ]

  rm(sc_exprs, bulk_dataset, matched_bulk)

  # Prepare X matrix
  x <- Matrix::t(Expression_bulk)

  # Step 3: Network construction with matrix optimizations
  cor.m <- if (p$network_class == 'bulk') {
    if (p$verbose) {
      ts_cli$cli_alert_info(
        "Constructing a gene-gene similarity by bulk data"
      )
    }
    stats::cor(x)
  } else {
    if (p$verbose) {
      ts_cli$cli_alert_info(
        "Constructing a gene-gene similarity by single cell data"
      )
    }
    # Use matrix operations for efficient correlation
    if (!inherits(Expression_cell, 'sparseMatrix')) {
      Expression_cell <- Matrix::Matrix(
        Expression_cell,
        sparse = TRUE
      )
    }
    scPAS::sparse.cor(Matrix::t(Expression_cell))
  }

  # Network construction
  cor.m[cor.m < 0] <- 0
  SNN <- Seurat::FindNeighbors(
    1 - cor.m,
    distance.matrix = TRUE,
    verbose = p$verbose
  )
  Network <- as.matrix(SNN$snn)
  diag(Network) <- 0
  Network <- (Network > 0.2) * 1 # binarization

  # Clean up
  rm(cor.m, SNN)
  gc(verbose = FALSE)

  # Step 4: Model optimization with purrr functional programming
  if (p$verbose) {
    ts_cli$cli_alert_info(
      "Optimizing the network-regularized sparse regression model"
    )
  }

  # Prepare Y based on family using purrr pattern matching
  y <- scPAS::prepare_phenotype(
    phenotype = phenotype,
    family = p$family,
    tag = p$label_type,
    verbose = p$verbose
  )

  model <- scPAS::optimize_model(
    x = x,
    y = y,
    Network = Network,
    alpha = alpha,
    cutoff = cutoff,
    family = p$family,
    seed = p$seed,
    verbose = p$verbose
  )
  # Step 5: Risk score calculation
  if (p$verbose) {
    ts_cli$cli_alert_info("Calculating quantified risk scores...")
  }

  # Sparse matrix scaling and risk calculation
  FastSparseRowScale <- utils::getFromNamespace("FastSparseRowScale", "Seurat")
  scaled_exp <- FastSparseRowScale(
    Expression_cell,
    display_progress = FALSE
  )
  scaled_exp[is.na(scaled_exp)] <- 0
  scaled_exp <- Matrix::Matrix(scaled_exp) # Probably a dgCMatrix
  # Matrix multiplication for risk scores
  risk_score <- Matrix::crossprod(scaled_exp, model$Coefs)

  # Step 6: Permutation test
  if (p$verbose) {
    ts_cli$cli_alert_info(
      "Qualitative identification by permutation test program with {.val {permutation_times}} times random perturbations..."
    )
  }

  bg_stats <- scPAS::perform_permutation_test(
    scaled_exp = scaled_exp,
    Coefs = model$Coefs,
    permutation_times = permutation_times,
    independent = independent,
    FDR.threshold = FDR_threshold,
    seed = p$seed
  )

  # Z-score calculation
  risk_score_df <- scPAS::calculate_risk_score(
    risk_score = risk_score,
    mean.background = bg_stats$mean.background,
    sd.background = bg_stats$sd.background,
    FDR.threshold = FDR_threshold,
    cell_names = colnames(Expression_cell)
  )

  sc_data <- SigBridgeRUtils::AddMisc(
    sc_data,
    scPAS = props(p$cache_config),
    cover = FALSE
  )

  sc_data$scPAS_RS <- risk_score_df$raw_score
  sc_data$scPAS_NRS <- risk_score_df$Z.statistics
  sc_data$scPAS_Pvalue <- risk_score_df$p.value
  sc_data$scPAS_FDR <- risk_score_df$FDR
  sc_data$scPAS <- risk_score_df$cell_label

  if (p$verbose) {
    ts_cli$cli_alert_success(
      cli::col_green("scPAS screening done.")
    )
  }

  list(
    scRNA_data = sc_data
  )
}
