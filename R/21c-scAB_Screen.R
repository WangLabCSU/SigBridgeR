#' @title Perform scAB Screening Analysis
#'
#' @description
#' Implements the scAB algorithm to identify phenotype-associated cell subpopulations
#' in single-cell RNA-seq data by integrating matched bulk expression and phenotype
#' information. Uses non-negative matrix factorization (NMF) with dual regularization
#' for phenotype association and cell-cell similarity.
#'
#' @inheritParams Screen
#' @param alpha Coefficient of phenotype regularization (default=`0.005`).
#' @param alpha_2 Coefficent of cell-cell similarity regularization (default=`0.005`).
#' @param k_max The maximum rank value to consider in the search. Must be at least 2. Defaults to `20`.
#' @param cross_k Number of folds for k-fold cross-validation. Defaults to `5`.
#' @param repeat_times The number of repeated NMF runs for each candidate rank to account for random initialization variability. Defaults to `10`.
#' @param maxiter Maximum number of iterations for NMF (default=2000).
#' @param tred Z-score threshold in finding subsets (default=2).
#' @param ... Additional arguments. Currently supports:
#'    - `verbose`: Logical indicating whether to print progress messages. Defaults to `TRUE`.
#'    - `seed`: For reproducibility, default is `123L`
#'    - `parallel`: Logical indicating whether to use parallel processing. Defaults to `FALSE`.
#'    - Other arguments are passed to `scAB::create_scAB.v5()`
#'    - `assay`: Character specifying the assay to use. Defaults to `"RNA"`.
#'    - `load_cache`: Cache directory path for loading precomputed scAB objects.
#'      Supports root-level, cache-level, or parent-level paths. See [CacheSetHere()].
#'    - `save_cache`: Cache directory path for saving scAB intermediate results.
#'      Supports root-level or parent-level paths. See [CacheSetHere()].
#'
#
#'
#' @return A list containing:
#' \describe{
#'   \item{scRNA_data}{Filtered Seurat object with selected cells}
#'   \item{scAB_result}{scAB screening result}
#' }
#'
#' @references
#' Zhang Q, Jin S, Zou X. scAB detects multiresolution cell states with clinical significance by integrating single-cell genomics and bulk sequencing data. Nucleic Acids Research. 2022 Nov 28;50(21):12112–30.
#'
#' @section LICENSE:
#' Licensed under the GNU General Public License version 3 (GPL-3.0).
#' A copy of the license is available at <https://www.gnu.org/licenses/gpl-3.0.en.html>.
#'
#' @examples
#' \dontrun{
#' # Binary phenotype example
#' result <- DoscAB(
#'   matched_bulk = bulk_matrix,
#'   sc_data = seurat_obj,
#'   phenotype = clinical_df,
#'   label_type = "disease_status",
#'   phenotype_class = "binary",
#'   alpha = 0.005,
#'   alpha_2 = 0.005,
#'   maxiter = 2000,
#'   tred = 2
#' )
#' }
#'
#' @export
#' @family screen_method
#' @family scAB
#'
DoscAB <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = "scAB",
  phenotype_class = c("binary", "survival"),
  alpha = c(0.005, NULL),
  alpha_2 = c(0.005, NULL),
  k_max = 20L,
  cross_k = 5L,
  repeat_times = 10L,
  maxiter = 2000L,
  tred = 2L,
  ...
) {
  check_installed("scAB", action = \(pkg, ...) {
    check_installed("pak")
    pak::pak("Exceret/scAB")
  })

  chk::chk_is(sc_data, "Seurat")
  chk::chk_character(label_type)
  phenotype_class <- arg_match(phenotype_class)
  chk::chk_number(maxiter)
  chk::chk_number(tred)
  # scAB can't tolerate NA
  chk::chk_not_any_na(matched_bulk)
  chk::chk_not_any_na(phenotype)

  # scAB is more strict than Scissor and scPAS
  if (phenotype_class == "survival") {
    if (!all(rownames(phenotype) == colnames(matched_bulk))) {
      cli::cli_abort(c(
        "x" = "Please check the rownames of {.var phenotype} and colnames of {.var bulk_dataset}, they should be the same."
      ))
    }
  } else {
    if (!all(names(phenotype) == colnames(matched_bulk))) {
      cli::cli_abort(c(
        "x" = "Please check the names of {.var phenotype} and colnames of {.var bulk_dataset}, they should be the same."
      ))
    }
  }

  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption("verbose")
  seed <- dots$seed %||% SigBridgeRUtils::getFuncOption("seed")
  assay <- dots$assay %||% "RNA"
  load_cache <- dots$load_cache
  save_cache <- dots$save_cache

  # -- build cache params ----------------------------------------------------
  cache_params <- ScreenMethodConfig(
    method_name = "scAB",
    phenotype_class = phenotype_class,
    label_type = label_type,
    param = fn_fmls() # NULL is OK; stores the parameters
  )

  # -- load mode: restore cached scAB_obj and k ------------------------------
  if (!is.null(load_cache)) {
    cache <- CacheSysCall(
      mode = "load",
      path = load_cache,
      cache = cache_params,
      verbose = verbose,
      timestamp = dots$timestamp,
    )

    scAB_obj <- cache$scAB_obj
    k <- cache$k

    if (verbose) {
      ts_cli$cli_alert_info(
        cli::col_green("Loaded scAB cache from {.path {cache_dir}}")
      )
    }
  } else {
    # -- normal flow: create scAB_obj and select k --------------------------
    if (verbose) {
      ts_cli$cli_alert_info(cli::col_green("Start scAB screening"))
    }

    scAB_obj <- scAB::create_scAB.v5(
      Object = sc_data,
      bulk_dataset = matched_bulk,
      phenotype = phenotype,
      method = phenotype_class,
      verbose = verbose,
      assay = assay
    )

    if (any(scAB_obj$X < 0)) {
      cli::cli_warn(
        "Found negative values in `X` after correlation, truncating to 0"
      )
      scAB_obj$X[scAB_obj$X < 0L] <- 0L
    }

    if (verbose) {
      ts_cli$cli_alert_info("Selecting K")
    }

    k <- scAB::select_K.optimized(
      Object = scAB_obj,
      K_max = k_max,
      repeat_times = repeat_times,
      maxiter = maxiter, # default in scAB
      seed = seed,
      verbose = verbose
    )

    # -- save mode: persist scAB_obj and k to cache --------------------------
    if (!is.null(save_cache)) {
      cache_dir <- CacheSysCall(
        mode = "save",
        path = save_cache,
        cache = cache_params,
        verbose = verbose,
        timestamp = dots$timestamp
      )
    }
  }

  if (verbose) {
    ts_cli$cli_alert_info(
      "Run NMF with phenotype and cell-cell similarity regularization at K = {.val {k}}"
    )
  }

  # Find optimal alpha and alpha_2
  if (
    is.null(alpha) ||
      is.null(alpha_2) ||
      length(alpha) > 1 ||
      length(alpha_2) > 1
  ) {
    para_list <- scAB::select_alpha.optimized(
      Object = scAB_obj,
      method = phenotype_class,
      K = k,
      cross_k = cross_k,
      para_1_list = alpha %||% c(0.01, 0.005, 0.001),
      para_2_list = alpha_2 %||% c(0.01, 0.005, 0.001),
      seed = seed,
      parallel = parallel,
      verbose = verbose
    )

    alpha <- para_list$para$alpha_1
    alpha_2 <- para_list$para$alpha_2
  }

  scAB_result <- scAB::scAB.optimized(
    Object = scAB_obj,
    K = k,
    alpha = alpha,
    alpha_2 = alpha_2,
    maxiter = maxiter,
    convergence_threshold = 1e-05
  )

  if (verbose) {
    ts_cli$cli_alert_info("Screening cells")
  }

  sc_data <- scAB::findSubset.optimized(
    Object = sc_data,
    scAB_Object = scAB_result,
    tred = tred
  )
  sc_data <- SigBridgeRUtils::AddMisc(
    seurat_obj = sc_data,
    scAB = cache_params,
    cover = FALSE
  )

  if (verbose) {
    ts_cli$cli_alert_info(
      cli::col_green("scAB screening done")
    )
  }

  list(scRNA_data = sc_data, scAB_result = scAB_result)
}
