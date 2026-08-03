#' @title Data-Driven Selection of Single-Cell Normalization Methods
#' @description
#' A quantitative framework for selecting optimal normalization strategies
#' (e.g., SCTransform vs. LogNormalization) based on diagnostic metrics rather
#' than heuristics. Evaluates three critical dimensions of preprocessing quality:
#' 1. **Variance stabilization**: Decoupling of mean-variance relationship
#'    in normalized expression (lower correlation = better).
#' 2. **Biological signal retention**: Preservation of known marker genes
#'    within highly variable genes (higher retention = better).
#' 3. **Dropout robustness**: Removal of technical dropout bias from
#'    normalized values (lower correlation with dropout rate = better).
#'
#' Methods are ranked using a weighted composite score. Designed for head-to-head
#' comparison of preprocessed Seurat objects.
#'
#' @param ... Named arguments where each value is a `Seurat` object representing
#'   a distinct preprocessing strategy (e.g., `SCT = sct_obj, Log = log_obj`).
#'   **Requirements**:
#'   * Must be named (names become method identifiers)
#'   * Must contain normalized data in the `data` slot of the specified assay/layer
#'   * Must have identical cell counts (for fair comparison)
#' @param subset_size Integer. Number of cells to subsample for diagnostics.
#'   If missing or empty, defaults to \code{min(10000, total_cells)}. Smaller subsets
#'   accelerate computation with minimal accuracy loss for large datasets (>50k cells).
#' @param known_hvgs Named list of canonical marker genes for cell types. Format:
#'   \code{list(cell_type1 = c("geneA", "geneB"), cell_type2 = c("geneC", ...))}.
#'   Used to compute \code{marker_retention} metric. If \code{NULL} (default), this
#'   metric is omitted from scoring.
#' @param n_hvgs Integer. Number of top highly variable genes to evaluate for marker
#'   retention. Default: \code{2000L}. Only relevant when \code{known_hvgs} is provided.
#' @param low_expressed_thresh Numeric (0–1). Quantile threshold for filtering lowly
#'   expressed genes during variance-mean correlation calculation. Genes below this
#'   quantile of mean expression are excluded. Default: \code{0.2} (bottom 20% excluded).
#' @param weight Named numeric vector specifying weights for composite scoring. Must
#'   contain exactly these components summing to 1:
#'
#'   * **`variance_stability`**: Weight for variance-mean decoupling (default: 0.4)
#'   * **`marker_signal`**: Weight for marker gene retention (default: 0.35)
#'   * **`dropout_robustness`**: Weight for dropout bias removal (default: 0.25)
#'
#' @return A list containing:
#'
#' * **`metrics`**: A `data.table` with diagnostic metrics per method:
#'   * `variance_mean_cor`: Pearson correlation between log10(mean) and
#'     log10(variance) of normalized expression (lower = better)
#'   * `marker_retention`: Proportion of known markers in top HVGs
#'     (higher = better; NA if `known_hvgs` not provided)
#'   * `mean_dropout_residual`: Absolute Spearman correlation between
#'     dropout rate (from counts) and normalized means (lower = better)
#'   * `composite_score`: Weighted combination of normalized metrics (0–1 scale)
#'   * `rank`: Method ranking (1 = best)
#' * **`recommendation`**: Character string naming the top-ranked method
#' * **`plots`**: (Optional) Diagnostic visualizations if `ggplot2` available
#'
#' @section Workflow:
#' 1. Validates input objects (naming, cell count consistency, data slot presence)
#' 2. Subsamples cells (if needed) for computational efficiency
#' 3. Computes three core metrics per method:
#'    * **Variance stabilization**: Correlation between log-transformed mean and
#'      variance of normalized expression
#'    * **Marker retention**: Overlap between user-provided markers and top HVGs
#'    * **Dropout robustness**: Correlation between gene dropout rate (from counts)
#'      and normalized expression means
#' 4. Normalizes metrics to 0-1 scale (inverting where lower=better)
#' 5. Computes weighted composite score and ranks methods
#' 6. Returns top recommendation with full diagnostic report
#'
#' @section Notes:
#' * Requires *pre-normalized* Seurat objects—this function **does not**
#'   perform normalization itself. Users must first run `SCTransform()`,
#'   `NormalizeData()`, etc., and store results in the `data` slot.
#' * Marker retention metric is only computed when `known_hvgs` is provided.
#'   Without it, scoring relies solely on variance stabilization and dropout robustness.
#' * Weights should sum to 1 (not enforced but recommended for interpretable scores).
#'
#' @examples
#' \dontrun{
#' # Compare two preprocessed Seurat objects
#' sct_obj <- SCTransform(seurat_raw, verbose = FALSE)
#' log_obj <- NormalizeData(seurat_raw, normalization.method = "LogNormalize")
#'
#' result <- ChooseNormalization(
#'   SCT = sct_obj,
#'   Log = log_obj,
#'   subset_size = 5000,
#'   verbose = TRUE
#' )
#'
#' # Top recommendation
#' result$recommendation
#'
#' # Full metrics table
#' result$metrics[, .(method, composite_score, rank)]
#'
#' # With marker genes for refined scoring
#' markers <- list(
#'   T_cell = c("CD3D", "CD3E", "CD8A"),
#'   B_cell = c("CD79A", "MS4A1", "CD19"),
#'   Myeloid = c("CD14", "LYZ", "FCGR3A")
#' )
#'
#' result <- ChooseNormalization(
#'   SCT = sct_obj,
#'   Log = log_obj,
#'   known_hvgs = markers,
#'   weight = c(variance_stability = 0.3, marker_signal = 0.5, dropout_robustness = 0.2)
#' )
#' }
#' @export
ChooseNormalization <- function(
  ...,
  subset_size = integer(),
  known_hvgs = list(),
  n_hvgs = 2000L,
  low_expressed_thresh = 0.2,
  weight = c(
    variance_stability = 0.4,
    marker_signal = 0.35,
    dropout_robustness = 0.25
  )
) {
  dots <- list2(...)
  is_seurat <- vapply(
    X = dots,
    FUN = \(x) inherits(x, "Seurat"),
    FUN.VALUE = logical(1)
  )
  method_objects <- dots[is_seurat]
  method_names <- names(method_objects)[is_seurat]
  n_cells <- purrr::map_int(method_objects, ~ ncol(.x))

  assay <- dots$assay
  layer <- dots$layer

  ChooseNormalizationCheck(
    method_objects = method_objects,
    method_names = method_names,
    n_cells = n_cells,
    assay = assay,
    layer = layer,
    low_expressed_thresh = low_expressed_thresh,
    weight = weight
  )

  verbose <- dots$verbose %||% getFuncOption("verbose")
  seed <- dots$seed %||% getFuncOption("seed")

  set.seed(seed)

  # * subset seurat objects
  if (is.null(subset_size) || length(subset_size) == 0) {
    subset_size <- min(n_cells, 10000)

    if (verbose) {
      cli::cli_alert_info("Using {subset_size} cells")
    }
  }
  method_objects <- purrr::map(
    method_objects,
    ~ subset(.x, cells = sample(colnames(.x), subset_size))
  )

  if (verbose) {
    cli::cli_alert_info(
      "Comparing {length(method_objects)} method{?s}: {.pkg {method_names}}"
    )
  }

  # * ===== 2. Extract metrics for each method =====
  metrics <- purrr::imap(
    method_objects,
    ~ exec(
      ExtractMetrics,
      obj = .x,
      assay = assay,
      low_expressed_thresh = low_expressed_thresh,
      method_name = .y,
      ground_truth_markers = known_hvgs,
      n_hvgs = n_hvgs,
      SigBridgeRUtils::FilterArgs4Func(dots, Seurat::VariableFeatures)
    ),
    .progress = if (verbose) "Calculating" else FALSE
  ) |>
    data.table::rbindlist()

  # * ===== 3. Normalize metrics for composite scoring =====
  # Lower is better: variance_mean_cor, mean_dropout_residual
  # Higher is better: marker_retention, hvgs_with_signal
  composite_score <- variance_stability <- marker_signal <- dropout_robustness <- NULL # ease checking NOTE
  metrics[,
    c("variance_stability", "dropout_robustness") := lapply(
      .SD,
      normalize_metric,
      invert = TRUE
    ),
    .SDcols = c("variance_mean_cor", "mean_dropout_residual")
  ]
  metrics[,
    c("marker_signal") := lapply(.SD, normalize_metric),
    .SDcols = c("marker_retention")
  ]

  # Composite score (weighted average)
  metrics[,
    composite_score := variance_stability *
      (weight["variance_stability"]) +
      marker_signal * (weight["marker_signal"]) +
      dropout_robustness * (weight["dropout_robustness"])
  ]

  data.table::setorder(metrics, -composite_score)
  metrics[, rank := 1:.N] # Rank methods

  plots <- ChooseNormalizationViz(metrics)

  # ===== 6. Report results =====
  if (verbose) {
    best_method <- metrics$method[1]

    cli::cli_h3("Method Ranking (Composite Score)")
    cli::cli_text(paste0(
      "Top method: {.field {best_method}} "
    ))

    purrr::walk(1:min(5, nrow(metrics)), function(i) {
      cli::cli_text(
        "{.strong [{i}]} {metrics$method[i]}: {.val {round(metrics$composite_score[i], 3)}}"
      )
    })
    # Highlight key differentiators
    best_metrics <- metrics[method == best_method]
    runner_up <- if (nrow(metrics) > 1) metrics$method[2] else NULL

    if (!is.null(runner_up)) {
      diff_vm <- abs(
        best_metrics$variance_mean_cor -
          metrics$variance_mean_cor[metrics$method == runner_up]
      )
      if (diff_vm > 0.15) {
        cli::cli_text(
          "{cli::symbol$bullet} Variance stabilization: {.val {best_method}} shows substantially better mean-variance decoupling"
        )
      }

      if (
        !is.na(best_metrics$marker_retention) &&
          best_metrics$marker_retention > 0.15
      ) {
        cli::cli_text(
          "{cli::symbol$bullet} Biological signal: {.val {best_method}} retains more marker genes in HVGs"
        )
      }
    }
  }

  list(metrics = metrics, recommendation = best_method)
}

#' @keywords internal
ChooseNormalizationCheck <- function(
  method_objects,
  method_names,
  n_cells,
  assay = NULL,
  layer = NULL,
  low_expressed_thresh = 0.2,
  weight = c(
    variance_stability = 0.4,
    marker_signal = 0.35,
    dropout_robustness = 0.25
  ),
  ...
) {
  # more than 1 method provided
  chk::chk_length(x = method_objects, length = 1L, upper = 30L)
  # methods are named
  purrr::walk(
    method_names,
    ~ if (!nzchar(.x, keepNA = TRUE)) {
      Abort(
        "Options must be named",
        "e.g., SCT = sct_obj, Log = log_obj"
      )
    }
  )
  # all Seurat objects
  purrr::walk(
    method_objects,
    ~ chk::chk_is(x = .x, class = "Seurat")
  )
  # contains data slot
  purrr::walk(
    method_objects,
    ~ if (
      "data" %chin%
        methods::slotNames(SeuratObject::LayerData(
          object = .x,
          assay = assay,
          layer = layer
        ))
    ) {
      if (is.null(layer)) {
        Abort(
          "All objects must contain normalized data in {.cls {assay}} assay {.cls data} slot"
        )
      }

      Abort(
        "All objects must contain normalized data in {.cls {assay}} assay {.cls {layer}} {.cls data} slot"
      )
    }
  )
  # same cell number
  if (length(unique(n_cells)) > 1) {
    Abort(
      "Seurat objects contain different cell counts. Ensure comparable subsets",
      tips = "Detected: {n_cells}"
    )
  }
  chk::chk_range(low_expressed_thresh)
  chk::chk_named(
    weight,
    x_name = c("variance_stability", "marker_signal", "dropout_robustness")
  )
  chk::chk_vector(weight)
  chk::chk_length(weight, 3)
  if (sum(weight) != 1) {
    Abort("weight must sum to 1")
  }
}

#' @keywords internal
ExtractMetrics <- function(
  obj,
  assay = "RNA",
  low_expressed_thresh = 0.2,
  method_name = character(),
  ground_truth_markers = NULL,
  n_hvgs = 2000,
  ...
) {
  dots <- list2(...)
  # 1. Variance-mean correlation on normalized data
  assay_data <- SeuratObject::LayerData(
    object = obj,
    assay = assay,
    layer = "data"
  )
  means <- SigBridgeRUtils::rowMeans3(assay_data)
  vars <- SigBridgeRUtils::rowVars3(assay_data)

  # Filter lowly expressed genes (bottom 20%)
  valid <- means > stats::quantile(means, low_expressed_thresh, na.rm = TRUE) &
    vars > 0
  var_mean_cor <- if (sum(valid) > 100) {
    if (is_installed(c("cheapr", "WGCNA"))) {
      WGCNA::cor(
        cheapr::log10_(means[valid] + 1),
        cheapr::log10_(vars[valid] + 1e-6),
        use = "complete.obs"
      )
    } else {
      stats::cor(
        log10(means[valid] + 1),
        log10(vars[valid] + 1e-6),
        use = "complete.obs"
      )
    }
  } else {
    NA
  }

  # 2. Marker gene retention in HVGs
  if (!is.null(ground_truth_markers)) {
    all_markers <- unique(unlist(ground_truth_markers))
    hvgs_detected <- exec(
      Seurat::VariableFeatures,
      object = obj,
      !!!SigBridgeRUtils::FilterArgs4Func(dots, Seurat::VariableFeatures)
    )
    n_hvgs_detected <- length(hvgs_detected)
    hvgs <- utils::head(hvgs_detected, n_hvgs)
    if (length(hvgs) > 0) {
      marker_retention <- length(intersect(hvgs, all_markers)) / length(hvgs)
    }
  } else {
    marker_retention <- NA
    n_hvgs_detected <- length(exec(
      Seurat::VariableFeatures,
      object = obj,
      !!!SigBridgeRUtils::FilterArgs4Func(dots, Seurat::VariableFeatures)
    ))
  }

  # 3. Dropout residual analysis
  # Assess whether normalized values still correlate with detection rate
  counts <- SeuratObject::LayerData(
    object = obj,
    assay = assay,
    layer = "counts"
  )
  dropout_rate <- SigBridgeRUtils::rowMeans3(counts == 0)

  # Correlation between dropout rate and normalized expression
  # Lower = better (successful normalization removes technical dropout bias)
  dropout_residual <- if (!is_installed("WGCNA")) {
    stats::cor(
      dropout_rate,
      means,
      use = "complete.obs",
      method = "spearman"
    )
  } else {
    WGCNA::cor(
      dropout_rate,
      means,
      use = "complete.obs",
      method = "spearman"
    )
  }

  list(
    method = method_name,
    variance_mean_cor = var_mean_cor,
    marker_retention = marker_retention,
    mean_dropout_residual = abs(dropout_residual), # absolute value for scoring
    n_cells = ncol(obj),
    n_genes = nrow(obj),
    n_hvgs_detected = n_hvgs_detected
  )
}


#' Min-max normalizati
#' @keywords internal
normalize_metric <- function(x, invert = FALSE) {
  sd1 <- if (is_installed("cheapr")) {
    cheapr::sd_(x, na.rm = TRUE)
  } else {
    stats::sd(x, na.rm = TRUE)
  }

  if (sd1 < .Machine$double.eps) {
    return(rep(0.5, length(x)))
  }
  norm_x <- (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  if (invert) 1 - norm_x else norm_x
}

#' @keywords internal
ChooseNormalizationViz <- function(metrics_df) {
  check_installed(c("ggplot2", "patchwork"))

  method <- variance_mean_cor <- composite_score <- NULL # ease checking NOTE
  # 1. Variance-mean correlation comparison
  vm_plot <- ggplot2::ggplot(
    metrics_df,
    ggplot2::aes(
      x = method,
      y = variance_mean_cor,
      fill = method
    )
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_hline(
      yintercept = 0.3,
      linetype = "dashed",
      color = "gray50",
      linewidth = 0.5
    ) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(
      title = "Variance-Mean Correlation (Lower = Better)",
      subtitle = "Measures residual technical noise after normalization",
      y = "Pearson correlation",
      x = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, size = 10)
    )

  # 2. Composite score radar (simplified bar version)
  score_plot <- ggplot2::ggplot(
    metrics_df,
    ggplot2::aes(
      x = method,
      y = composite_score,
      fill = method
    )
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    ggplot2::labs(
      title = "Composite Quality Score (Higher = Better)",
      y = "Normalized score (0-1)",
      x = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, size = 10)
    )

  list(
    var_mean_plot = vm_plot,
    composite_score_plot = score_plot,
    combined_plot = patchwork::wrap_plots(vm_plot, score_plot, ncol = 2)
  )
}
