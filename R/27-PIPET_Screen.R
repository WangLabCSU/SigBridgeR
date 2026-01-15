#' @title Perform PIPET Screening Analysis
#' @description
#' Predicts cell subpopulations in single-cell data by matching expression profiles
#' to predefined marker gene templates using various distance/similarity metrics.
#' This function implements a template-based classification approach with permutation
#' testing for significance assessment.
#'
#'
#' @param matched_bulk Normalized bulk expression matrix (features × samples).
#'        Column names must match `phenotype` identifiers.
#' @param sc_data Seurat object containing single-cell RNA-seq data.
#' @param phenotype Clinical outcome data. Can be:
#'        - Vector: named with sample IDs
#'        - Data frame: with row names matching bulk columns
#' @param phenotype_class Analysis mode:
#'        - `"binary"`: Case-control design (e.g., responder/non-responder)
#'        - `"continuous"`: Continuous outcome (e.g.,)
#' @param group A character, name of one metadata column to group cells by (for example, orig.ident). The default value is `NULL`. In this case, screening will be performed on each group separately.
#' @param discretize_method \code{c("median", "kmeans", "custom")}. Discretization
#'   strategy for continuous phenotypes. Note: `"median"` is mapped internally to
#'   `"quantile"` (2-group quantile split). Default: `"median"`.
#' @param cutoff Numeric vector of length `n_group - 1`. Required only when
#'   \code{discretize_method = "custom"}. Defines interior breakpoints on the
#'   *normalized, log2-transformed scale* (i.e., after `scale(log2(x + 1))`).
#'   Must be sorted in ascending order.
#' @param label_type Character specifying phenotype label type (e.g., "SBS1", "time"), stored in `scRNA_data@misc`
#' @param log2FC In the DESeq differential expression analysis results, the cutoff value of log2FC. The default value is `1L`.
#' @param p_adjust In the DESeq differential expression analysis results, the cutoff value of adjust P. The default value is `0.05`.
#' @param show_log2FC Select whether to show log2 fold changes. The default value is `TRUE`.
#' @param freq_counts An integer, keep genes expressed in more than a certain number of cells. The default value is `NULL`, which means no filtering.
#' @param normalize Select whether to perform normalization of count data. The default value is `TRUE`.
#' @param scale Select whether to scale and center features in the dataset. The default value is `TRUE`.
#' @param nPerm An integer, number of permutations to do. The default value is `1000L`.
#' @param distance A character, the distance algorithm must be included in "cosine", "pearson", "spearman", "kendall","euclidean","maximum". default value is `NULL`, which means `"cosine"`.
#' @param ... Additional arguments to be passed to `PIPET.optimized`.
#' - seed: Random seed for reproducibility
#' - verbose: Whether to show progress messages
#' - parallel: Whether to use parallel processing, default is `FALSE`. future::plan() must be set before calling this function.
#'
#'
#' @family screen_method
#' @family PIPET
#' @export
DoPIPET <- function(
  matched_bulk,
  sc_data,
  phenotype,
  phenotype_class = c("binary", "continuous"),
  group = NULL,
  discretize_method = c("kmeans", "median", "custom"),
  cutoff = NULL,
  label_type = "PIPET",
  log2FC = 1L,
  p_adjust = 0.05,
  show_log2FC = TRUE,
  freq_counts = NULL,
  normalize = TRUE,
  scale = TRUE,
  nPerm = 1000L,
  distance = c(
    "cosine",
    "pearson",
    "spearman",
    "kendall",
    "euclidean",
    "maximum"
  ),
  ...
) {
  # * Input validation
  chk::chk_is(matched_bulk, c("matrix", "data.frame"))
  chk::chk_is(sc_data, "Seurat")
  chk::chk_character(label_type)
  phenotype_class <- SigBridgeRUtils::MatchArg(
    phenotype_class,
    c("binary", "continuous"),
    NULL
  )
  purrr::walk(c(show_log2FC, rm_NA, normalize, scale), chk::chk_flag)
  purrr::walk(
    c(group, distance),
    ~ {
      if (!is.null(.x)) chk::chk_character(.x)
    }
  )
  chk::chk_integer(nPerm)
  purrr::walk(c(log2FC, p_adjust), chk::chk_double)

  # * Default params
  ## * general params
  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose")
  seed <- dots$seed %||% getFuncOption("seed")
  parallel <- (dots$parallel %||% FALSE) &
    !inherits(future::plan("list")[[1]], "sequential")
  ## * PIPET params
  distance <- SigBridgeRUtils::MatchArg(
    distance,
    c(
      "cosine",
      "pearson",
      "spearman",
      "kendall",
      "euclidean",
      "maximum"
    )
  )

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green("Starting PIPET screen"))
    ts_cli$cli_alert_info("Creating marker genes from bulk data...")
  }

  phenotype_df <- PIPET::AdaptPheno(
    phenotype = phenotype,
    phenotype_type = phenotype_class,
    discretize_method = discretize_method,
    cutoff = cutoff
  )

  markers <- PIPET::Create_Markers2(
    bulk_data = matched_bulk,
    colData = phenotype_df,
    class_col = "class",
    log2FC = log2FC,
    p.adjust = p_adjust,
    show_log2FC = show_log2FC,
    verbose = verbose,
    seed = seed
  )

  # Run PIPET core algorithm
  if (verbose) {
    ts_cli$cli_alert_info("Running PIPET correlation analysis...")
  }

  pipet_result <- PIPET::PIPET(
    sc_data = sc_data,
    markers = markers,
    group = group,
    freq_counts = freq_counts,
    normalize = normalize,
    scale = scale,
    nPerm = nPerm,
    distance = distance,
    verbose = verbose,
    seed = seed,
    parallel = parallel
  )

  # Add results to Seurat object if applicable

  sc_data <- Seurat::AddMetaData(
    sc_data,
    metadata = pipet_result
  ) %>%
    AddMisc(
      PIPET_type = label_type,
      PIPET_params = list(
        phenotype_class = phenotype_class,
        group = group,
        discretize_method = if (length(discretize_method) > 1) {
          discretize_method[1]
        } else {
          discretize_method
        },
        cutoff = cutoff,
        log2FC = log2FC,
        p_adjust = p_adjust,
        show_log2FC = show_log2FC,
        freq_counts = freq_counts,
        normalize = normalize,
        scale = scale,
        nPerm = nPerm,
        distance = if (length(distance) > 1) distance[1] else distance
      ),
      cover = FALSE
    )

  if (verbose) {
    ts_cli$cli_alert_success(cli::col_green("PIPET screening done."))
  }

  list(
    scRNA_data = sc_data,
    markers = markers
  )
}
