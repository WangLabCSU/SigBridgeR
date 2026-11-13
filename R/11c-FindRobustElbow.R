#' @title Automatically determine optimal PCA dimensions using multiple robust methods
#'
#' @description
#' This function combines multiple statistical approaches to automatically determine
#' the optimal number of principal components (PCs) for downstream single-cell analysis.
#' It integrates variance-based heuristics, elbow detection algorithms, and provides
#' comprehensive visualization for result validation.
#'
#' @param obj A Seurat object that has PCA computed (after `RunPCA`) or a numeric matrix
#'            of cell embeddings (cells x PCs).
#' @param verbose Logical, if TRUE outputs detailed method results and creates
#'                visualization plot. If `FALSE` returns only the final dimension.
#' @param ndims Integer, maximum number of dimensions to consider (default: `50L`).
#'   Values exceeding the number of PCs stored in the `"pca"` reduction are
#'   capped at the available number of PCs.
#'
#' @return Integer, the recommended number of PCA dimensions for downstream analysis
#'
#' @details
#' The numerical work is delegated to the C++ routine `find_robust_elbow_cpp()`,
#' which receives the PCA cell embeddings as an S4-aware matrix view created by
#' `beachmat::initializeCpp()`. Columns are streamed through `tatami` (see
#' `Rtatami.h`), so the embeddings are never coerced into a dense base `matrix`
#' and only the requested `ndims` columns are touched. For each of those PCs the
#' routine computes the standard deviation, the explained variance, the
#' percentage of explained variance, and the cumulative percentage, which are
#' then fed into the six heuristics described below.
#'
#' The recommended dimension is the second largest value among the following
#' candidates:
#'
#' * cumulative variance above 90\%;
#' * cumulative variance above 80\%;
#' * last PC explaining more variance than the mean;
#' * last PC explaining more than two standard deviations of the variance;
#' * maximum second derivative of the variance curve;
#' * largest distance between the variance curve and the chord joining the first
#'   and the last PC;
#' * the fixed anchors 10, 15, 20, 25, 30 and 35.
#'
#' The fixed anchors are independent of the data, so the raw recommendation can
#' exceed the number of PCs that were actually computed. The final value is
#' therefore capped at the number of available PCs (i.e. `ndims`).
#'
#' When `verbose = TRUE` and the session is interactive, the per-method
#' candidates are plotted on top of the variance curve by the internal helper
#' `plot_robust_elbow()`.
#'
#' @examples
#' \dontrun{
#' # After running PCA on Seurat object
#' pbmc <- RunPCA(pbmc, npcs = 50)
#' optimal_dims <- FindRobustElbow(pbmc, verbose = TRUE)
#' pbmc <- FindNeighbors(pbmc, dims = 1:optimal_dims)
#' }
#'
#' @family single_cell_preprocess
#' @export
#'
FindRobustElbow <- function(
  obj,
  verbose = SigBridgeRUtils::getFuncOption("verbose") %||% TRUE,
  ndims = 50L
) {
  # Input validation
  chk::chk_flag(verbose, x_name = "verbose")
  chk::chk_count(ndims, x_name = "ndims")
  chk::chk_gt(ndims, 0L, x_name = "ndims")

  # PCA cell embeddings: cells x PCs
  embeddings <- if (inherits(obj, "Seurat")) {
    if (!"pca" %chin% names(obj)) {
      Abort(
        "PCA has not been computed on this Seurat object.",
        "Please run {.code Seurat::RunPCA()} first.",
        type = "[DATA ERROR]"
      )
    }
    obj[["pca"]]@cell.embeddings
  } else {
    obj
  }
  available <- ncol(embeddings)

  if (ndims > available) {
    if (verbose) {
      cli::cli_warn(
        "Only {available} PC{?s} available; using {.code ndims = {available}}."
      )
    }
    ndims <- available
  }

  # All numerical work happens in C++
  res <- find_robust_elbow_cpp(beachmat::initializeCpp(embeddings), ndims)
  final_dims <- res$final_dims

  if (verbose) {
    cli::cli_h3(cli::col_green("Method Results"))
    cli::cli_alert_info(
      "Method 1A (Cumulative Variance > 90%): 1:{res$cumulative_90}"
    )
    cli::cli_alert_info(
      "Method 1B (Cumulative Variance > 80%): 1:{res$cumulative_80}"
    )
    cli::cli_alert_info(
      "Method 1C (Variance > Mean): 1:{res$above_mean}"
    )
    cli::cli_alert_info(
      "Method 1D (Variance > 2*SD): 1:{res$above_2sd}"
    )
    cli::cli_alert_info(
      "Method 2 (Second Derivative): 1:{res$second_deriv}"
    )
    cli::cli_alert_info(
      "Method 3 (Distance-based): 1:{res$distance}"
    )
    cli::cli_alert_success(
      "Final Recommended Dimensions:  1:{final_dims}"
    )

    # Create comprehensive visualization
    if (is_interactive()) {
      print(plot_robust_elbow(res))
    }

    # Additional summary
    cli::cli_h3(cli::col_green("Summary"))
    cli::cli_alert_info(
      "Cumulative variance at {.val {final_dims}} PCs: {.val {round(res$cumulative_variance[final_dims], 1)}}%"
    )
    cli::cli_alert_info(
      "Variance explained by PC{final_dims}: {.val {round(res$pct_variance[final_dims], 2)}}%"
    )

    # if (cumulative_variance[final_dims] < 80) {
    #     cli::cli_warn(
    #         "Cumulative variance < 80%. Consider increasing dimensions."
    #     )
    # }
  }

  min(final_dims, available)
}

#' @title Visualize the candidates proposed by robust elbow detection
#'
#' @description
#' Internal helper used by [FindRobustElbow()] to draw the variance curve of the
#' PCA dimensions, the per-method candidate cut-offs, and their labels.
#'
#' @param res A list as returned by `find_robust_elbow_cpp()`, containing at
#'   least `pct_variance`, `cumulative_variance`, `cumulative_90`,
#'   `cumulative_80`, `above_mean`, `above_2sd`, `second_deriv` and `distance`.
#'
#' @return A `ggplot2` object.
#'
#' @keywords internal
#'
plot_robust_elbow <- function(res) {
  PC <- Variance <- NULL # ease checking NOTE

  max_variance <- max(res$pct_variance)

  plot_data <- data.table::data.table(
    PC = seq_along(res$pct_variance),
    Variance = res$pct_variance,
    Cumulative = res$cumulative_variance
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = PC, y = Variance)
  ) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_line(alpha = 0.8) +
    ggplot2::geom_vline(
      xintercept = res$cumulative_90,
      color = "#0a9696ff",
      linetype = "dashed",
      alpha = 0.8,
      linewidth = 1.2
    ) +
    ggplot2::geom_vline(
      xintercept = res$cumulative_80,
      color = "#d699ffff",
      linetype = "dashed",
      alpha = 0.8,
      linewidth = 1.2
    ) +
    ggplot2::geom_vline(
      xintercept = res$above_mean,
      color = "#63118aff",
      linetype = "dashed",
      alpha = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = res$above_2sd,
      color = "#193db3ff",
      linetype = "dashed",
      alpha = 0.8,
      linewidth = 1.2
    ) +
    ggplot2::geom_vline(
      xintercept = res$second_deriv,
      color = "#15b915ff",
      linetype = "dashed",
      alpha = 0.8,
      linewidth = 1.2
    ) +
    ggplot2::geom_vline(
      xintercept = res$distance,
      color = "#ffac30ff",
      linetype = "dashed",
      alpha = 0.8,
      linewidth = 1.2
    ) +
    ggplot2::annotate(
      "text",
      x = res$cumulative_90,
      y = max_variance * 0.5,
      label = "Cumulative Variance > 90%",
      color = "#0a9696ff",
      hjust = -0.1,
      fontface = "bold"
    ) +
    ggplot2::annotate(
      "text",
      x = res$cumulative_80,
      y = max_variance * 0.4,
      label = "Cumulative Variance > 80%",
      color = "#d699ffff",
      hjust = -0.1,
      fontface = "bold"
    ) +
    ggplot2::annotate(
      "text",
      x = res$above_mean,
      y = max_variance * 0.9,
      label = "Variance > Mean",
      color = "#63118aff",
      hjust = -0.1,
      fontface = "bold"
    ) +
    ggplot2::annotate(
      "text",
      x = res$above_2sd,
      y = max_variance * 0.8,
      label = "Variance > 2*SD",
      color = "#193db3ff",
      hjust = -0.1,
      size = 4L,
      fontface = "bold"
    ) +
    ggplot2::annotate(
      "text",
      x = res$second_deriv,
      y = max_variance * 0.7,
      label = "Second Derivative",
      color = "#15b915ff",
      hjust = -0.1,
      size = 4L,
      fontface = "bold"
    ) +
    ggplot2::annotate(
      "text",
      x = res$distance,
      y = max_variance * 0.6,
      label = "Distance-based",
      color = "#ffac30ff",
      hjust = -0.1,
      size = 4L,
      fontface = "bold"
    ) +
    ggplot2::labs(
      title = "Robust Elbow Detection for PCA Dimension Selection",
      subtitle = "Multiple methods combined for optimal PC selection",
      x = "Principal Component",
      y = "Percentage of Variance Explained"
    ) +
    ggplot2::theme_minimal()
}
