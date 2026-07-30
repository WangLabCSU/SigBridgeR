# ---- 2. Do PIPET ----

#' Validate and Prepare DoPIPET Parameters
#'
#' @description
#' Internal helper that handles package installation checks, input validation,
#' and default-value resolution for [DoPIPET()].
#'
#' @param matched_bulk,sc_data,phenotype,phenotype_class,only_pos_marker,group
#'   Forwarded from [DoPIPET()].
#' @param discretize_method,cutoff,label_type,marker_finder,log2FC,p_adjust
#'   Forwarded from [DoPIPET()].
#' @param show_log2FC,freq_counts,normalize,scale,nPerm,distance
#'   Forwarded from [DoPIPET()].
#' @param ... Additional dots forwarded from [DoPIPET()].
#'
#' @return A named list with elements: `phenotype_class`, `marker_finder`,
#'   `distance`, `verbose`, `seed`, `parallel`, `assay`, `load_cache`,
#'   `save_cache`, `cache_config`, `dots`.
#'
#' @keywords internal
#' @family PIPET
ValidatePIPETParams <- function(
  matched_bulk,
  sc_data,
  phenotype,
  phenotype_class,
  only_pos_marker,
  group,
  discretize_method,
  cutoff,
  label_type,
  marker_finder,
  log2FC,
  p_adjust,
  show_log2FC,
  freq_counts,
  normalize,
  scale,
  nPerm,
  distance,
  ...
) {
  # -- package checks -------------------------------------------------------
  check_installed("PIPET", action = \(pkg, ...) {
    check_installed("pak")
    pak::pak("Exceret/PIPET")
  })

  # -- input validation -----------------------------------------------------
  chk::chk_is(matched_bulk, c("matrix", "data.frame"))
  chk::chk_is(sc_data, "Seurat")
  chk::chk_character(label_type)
  phenotype_class <- SigBridgeRUtils::MatchArg(
    phenotype_class,
    c("binary", "continuous", "survival"),
    NULL
  )
  purrr::walk(c(show_log2FC, normalize, scale), chk::chk_flag)
  purrr::walk(
    c(group),
    ~ if (!is.null(.x)) chk::chk_character(.x)
  )
  chk::chk_integer(c(nPerm, freq_counts))
  purrr::walk(c(log2FC, p_adjust), chk::chk_double)
  marker_finder <- SigBridgeRUtils::MatchArg(
    marker_finder,
    c("limma", "DESeq2")
  )

  distance <- SigBridgeRUtils::MatchArg(
    distance,
    c(
      "cosine",
      "pearson",
      "spearman",
      "kendall",
      "euclidean",
      "maximum"
    ),
    NULL
  )

  # -- process dots ---------------------------------------------------------
  dots <- list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose")
  seed <- dots$seed %||% getFuncOption("seed")
  parallel <- (dots$parallel %||% FALSE)
  assay <- dots$assay %||% "RNA"
  load_cache <- dots$load_cache
  save_cache <- dots$save_cache

  # -- build cache config ---------------------------------------------------
  cache_config <- ScreenMethodConfig(
    method_name = "PIPET",
    param = get_env_vars(exclude = c("matched_bulk", "sc_data", "phenotype")),
    phenotype_class = phenotype_class,
    label_type = label_type
  )

  get_env_vars(exclude = c("matched_bulk", "sc_data", "phenotype"))
}

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
#'        - `"continuous"`: Continuous outcome (e.g., age, size)
#'        - `"survival"`: Patient survival
#' @param only_pos_marker A logical value. If `TRUE`, only upregulated marker genes
#'        (those with positive log2 fold change) will be retained. If `FALSE` (default),
#'        both upregulated and downregulated markers are kept.
#' @param group A character, name of one metadata column to group cells by (for example, orig.ident). The default value is `NULL`. In this case, screening will be performed on each group separately.
#' @param discretize_method \code{c("median", "kmeans", "custom")}. Discretization
#'   strategy for continuous phenotypes. Note: `"median"` is mapped internally to
#'   `"quantile"` (2-group quantile split). Default: `"kmeans"`.
#' @param cutoff Numeric vector of length `n_group - 1`. Required only when
#'   \code{discretize_method = "custom"}. Defines interior breakpoints on the
#'   *normalized, log2-transformed scale* (i.e., after `scale(log2(x + 1))`).
#'   Must be sorted in ascending order.
#' @param label_type Character specifying phenotype label type (e.g., "SBS1", "time"), stored in `scRNA_data@misc`
#' @param marker_finder A character, the marker finder method. The default value is `"limma"`.
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
#' - assay: The assay to use, default is `"RNA"`
#'
#'
#' @family screen_method
#' @family PIPET
#' @export
DoPIPET <- function(
  matched_bulk,
  sc_data,
  phenotype,
  phenotype_class = c("binary", "continuous", "survival"),
  only_pos_marker = FALSE,
  group = NULL,
  discretize_method = c("kmeans", "median", "custom"),
  cutoff = NULL,
  label_type = "PIPET",
  marker_finder = c("limma", "DESeq2"),
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
  # -- validate & prepare all parameters -----------------------------------
  p <- exec(ValidatePIPETParams, !!!fn_fmls())

  if (p$verbose) {
    ts_cli$cli_alert_info(cli::col_green("Starting PIPET screen"))
  }

  # -- load mode: restore cached markers ------------------------------------
  if (!is.null(p$load_cache)) {
    cache <- CacheSysCall(
      mode = "load",
      path = p$load_cache,
      cache = p$cache_config,
      verbose = p$verbose,
      timestamp = p$dots$timestamp,
    )

    markers <- cache$markers
    rm(cache)

    if (p$verbose) {
      ts_cli$cli_alert_info(
        cli::col_green("Loaded PIPET markers from cache")
      )
    }
  } else {
    # -- normal flow: create markers ----------------------------------------
    if (p$verbose) {
      ts_cli$cli_alert_info("Creating markers from bulk data")
    }

    phenotype_df <- PIPET::AdaptPheno(
      phenotype = phenotype,
      phenotype_type = p$phenotype_class,
      discretize_method = discretize_method,
      cutoff = cutoff
    )

    # a data.frame
    markers <- if (p$marker_finder == "limma") {
      PIPET::Create_Markers2(
        bulk_data = matched_bulk,
        colData = phenotype_df,
        class_col = "class",
        lg2FC = log2FC,
        p.adjust = p_adjust,
        show_log2FC = show_log2FC,
        verbose = p$verbose,
        seed = p$seed
      )
    } else {
      PIPET::Create_Markers(
        bulk_data = matched_bulk,
        colData = phenotype_df,
        class_col = "class",
        lg2FC = log2FC,
        p.adjust = p_adjust,
        show_log2FC = show_log2FC,
        verbose = p$verbose,
        seed = p$seed
      )
    }

    if (only_pos_marker) {
      markers <- markers[markers$log2FoldChange > 0, ]
    }

    # -- save mode: persist markers to cache --------------------------------
    if (!is.null(p$save_cache)) {
      cache <- ScreenMethodCache(
        cache_path = p$save_cache,
        cache_config_path = file.path(p$save_cache, "cache_config.json"),
        cache_data = list(markers = markers),
        screen_method_config = p$cache_config
      )
      CacheSysCall(
        mode = "save",
        path = p$save_cache,
        cache = cache,
        verbose = p$verbose,
        timestamp = p$dots$timestamp
      )
    }
  }

  # Run PIPET core algorithm
  if (p$verbose) {
    ts_cli$cli_alert_info("Running PIPET correlation analysis")
  }

  pipet_result <- PIPET::PIPET(
    sc_data = sc_data,
    markers = markers,
    group = group,
    freq_counts = freq_counts,
    normalize = normalize,
    scale = scale,
    nPerm = nPerm,
    distance = p$distance,
    verbose = p$verbose,
    seed = p$seed,
    parallel = p$parallel,
    assay = p$assay
  )

  if (is.null(pipet_result)) {
    cli::cli_abort(c(
      "x" = "PIPET screening failed.",
      ">" = "Try different parameters"
    ))
  }

  # Add results to Seurat object
  sc_data <- Seurat::AddMetaData(
    sc_data,
    metadata = pipet_result
  )
  sc_data <- AddMisc(
    seurat_obj = sc_data,
    PIPET = props(p$cache_config),
    cover = FALSE
  )

  if (p$verbose) {
    ts_cli$cli_alert_success(cli::col_green("PIPET screening done."))
  }

  list(
    scRNA_data = sc_data,
    markers = markers
  )
}
