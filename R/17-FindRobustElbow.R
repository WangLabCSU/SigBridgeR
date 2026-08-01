#' @title Automatically determine optimal PCA dimensions using multiple robust methods
#'
#' @description
#' This function combines multiple statistical approaches to automatically determine
#' the optimal number of principal components (PCs) for downstream single-cell analysis.
#' It integrates variance-based heuristics, elbow detection algorithms, and provides
#' comprehensive visualization for result validation.
#'
#' @param obj A Seurat object that has PCA computed (after `RunPCA`)
#' @param verbose Logical, if TRUE outputs detailed method results and creates
#'                visualization plot. If `FALSE` returns only the final dimension.
#' @param ndims Integer, maximum number of dimensions to consider (default: `50L`)
#'
#' @return Integer, the recommended number of PCA dimensions for downstream analysis
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
  if (!"pca" %chin% names(obj)) {
    Abort(
      "PCA has not been computed on this Seurat object.",
      "Please run {.code RunPCA()} first."
    )
  }

  # Extract PCA standard deviations
  stdev <- obj[["pca"]]@stdev[seq_len(ndims)]

  # Calculate variance metrics
  variance <- stdev^2
  pct_variance <- variance / sum(variance) * 100
  cumulative_variance <- cumsum(pct_variance)

  # Method 1: Variance-based heuristics
  FindElbowVariance <- function(cumulative_variance, pct_variance) {
    method1_results <- list()

    # 1A: Cumulative variance > 90%
    dims_90pct <- which(cumulative_variance > 90)[1]
    method1_results$cumulative_90 <- ifelse(
      is.na(dims_90pct),
      ndims,
      dims_90pct
    )

    # 1B: Cumulative variance > 80%
    dims_80pct <- which(cumulative_variance > 80)[1]
    method1_results$cumulative_80 <- ifelse(
      is.na(dims_80pct),
      ndims,
      dims_80pct
    )

    # 1C: PCs explaining more than mean variance
    dims_above_mean <- which(pct_variance > mean(pct_variance))
    method1_results$above_mean <- ifelse(
      length(dims_above_mean) > 0,
      max(dims_above_mean),
      10
    )

    # 1D: PCs explaining more than 2*SD of variance
    threshold_2sd <- 2 * stats::sd(pct_variance)
    dims_above_2sd <- which(pct_variance > threshold_2sd)
    method1_results$above_2sd <- ifelse(
      length(dims_above_2sd) > 0,
      max(dims_above_2sd),
      10
    )

    return(method1_results)
  }
  method1_results <- FindElbowVariance(cumulative_variance, pct_variance)

  # Method 2: Second derivative elbow detection
  FindElbowDeriv <- function(variance, ndims) {
    first_deriv <- diff(variance)
    second_deriv <- diff(first_deriv)
    elbow_point <- which.max(abs(second_deriv)) + 1
    return(min(elbow_point, ndims))
  }

  method2_final <- FindElbowDeriv(variance, ndims)

  # Method 3: Robust distance-based elbow detection
  FindElbowDist <- function(variance, ndims) {
    # Calculate distances from each point to the line connecting first and last points
    line_start <- c(1, variance[1])
    line_end <- c(ndims, variance[ndims])

    distances <- vapply(
      seq_len(ndims),
      function(i) {
        point <- c(i, variance[i])
        # Calculate perpendicular distance from point to line
        numerator <- abs(
          (line_end[2] - line_start[2]) *
            point[1] -
            (line_end[1] - line_start[1]) * point[2] +
            line_end[1] * line_start[2] -
            line_end[2] * line_start[1]
        )
        denominator <- sqrt(
          (line_end[2] - line_start[2])^2 +
            (line_end[1] - line_start[1])^2
        )
        numerator / denominator
      },
      FUN.VALUE = numeric(1)
    )

    elbow_point <- which.max(distances)
    return(elbow_point)
  }

  method3_final <- FindElbowDist(variance, ndims)

  all_methods <- c(
    unlist(method1_results),
    method2_final,
    method3_final,
    10,
    15,
    20,
    25,
    30,
    35
  )
  # Second largest value is the final recommended dimension
  final_dims <- sort(all_methods, decreasing = TRUE)[2]

  if (verbose) {
    cli::cli_h3(cli::col_green("Method Results"))
    cli::cli_alert_info(
      "Method 1A (Cumulative Variance > 90%): 1:{method1_results$cumulative_90}"
    )
    cli::cli_alert_info(
      "Method 1B (Cumulative Variance > 80%): 1:{method1_results$cumulative_80}"
    )
    cli::cli_alert_info(
      "Method 1C (Variance > Mean): 1:{method1_results$above_mean}"
    )
    cli::cli_alert_info(
      "Method 1D (Variance > 2*SD): 1:{method1_results$above_2sd}"
    )
    cli::cli_alert_info(
      "Method 2 (Second Derivative): 1:{method2_final}"
    )
    cli::cli_alert_info(
      "Method 3 (Distance-based): 1:{method3_final}"
    )
    cli::cli_alert_success(
      "Final Recommended Dimensions:  1:{final_dims}"
    )

    # Create comprehensive visualization
    if (is_interactive()) {
      check_installed("ggplot2")
      PC <- Variance <- NULL # ease checking NOTE

      plot_data <- data.table::data.table(
        PC = seq_len(ndims),
        Variance = pct_variance,
        Cumulative = cumulative_variance
      )

      p <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = PC, y = Variance)
      ) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::geom_line(alpha = 0.8) +
        ggplot2::geom_vline(
          xintercept = method1_results$cumulative_90,
          color = "#0a9696ff",
          linetype = "dashed",
          alpha = 0.8,
          linewidth = 1.2
        ) +
        ggplot2::geom_vline(
          xintercept = method1_results$cumulative_80,
          color = "#d699ffff",
          linetype = "dashed",
          alpha = 0.8,
          linewidth = 1.2
        ) +
        ggplot2::geom_vline(
          xintercept = method1_results$above_mean,
          color = "#63118aff",
          linetype = "dashed",
          alpha = 0.8
        ) +
        ggplot2::geom_vline(
          xintercept = method1_results$above_2sd,
          color = "#193db3ff",
          linetype = "dashed",
          alpha = 0.8,
          linewidth = 1.2
        ) +
        ggplot2::geom_vline(
          xintercept = method2_final,
          color = "#15b915ff",
          linetype = "dashed",
          alpha = 0.8,
          linewidth = 1.2
        ) +
        ggplot2::geom_vline(
          xintercept = method3_final,
          color = "#ffac30ff",
          linetype = "dashed",
          alpha = 0.8,
          linewidth = 1.2
        ) +
        ggplot2::annotate(
          "text",
          x = method1_results$cumulative_90,
          y = max(pct_variance) * 0.5,
          label = "Cumulative Variance > 90%",
          color = "#0a9696ff",
          hjust = -0.1,
          fontface = "bold"
        ) +
        ggplot2::annotate(
          "text",
          x = method1_results$cumulative_80,
          y = max(pct_variance) * 0.4,
          label = "Cumulative Variance > 80%",
          color = "#d699ffff",
          hjust = -0.1,
          fontface = "bold"
        ) +
        ggplot2::annotate(
          "text",
          x = method1_results$above_mean,
          y = max(pct_variance) * 0.9,
          label = "Variance > Mean",
          color = "#63118aff",
          hjust = -0.1,
          fontface = "bold"
        ) +
        ggplot2::annotate(
          "text",
          x = method1_results$above_2sd,
          y = max(pct_variance) * 0.8,
          label = "Variance > 2*SD",
          color = "#193db3ff",
          hjust = -0.1,
          size = 4,
          fontface = "bold"
        ) +
        ggplot2::annotate(
          "text",
          x = method2_final,
          y = max(pct_variance) * 0.7,
          label = "Second Derivative",
          color = "#15b915ff",
          hjust = -0.1,
          size = 4,
          fontface = "bold"
        ) +
        ggplot2::annotate(
          "text",
          x = method3_final,
          y = max(pct_variance) * 0.6,
          label = "Distance-based",
          color = "#ffac30ff",
          hjust = -0.1,
          size = 4,
          fontface = "bold"
        ) +
        ggplot2::labs(
          title = "Robust Elbow Detection for PCA Dimension Selection",
          subtitle = "Multiple methods combined for optimal PC selection",
          x = "Principal Component",
          y = "Percentage of Variance Explained"
        ) +
        ggplot2::theme_minimal()

      methods::show(p)
    }

    # Additional summary
    cli::cli_h3(cli::col_green("Summary"))
    cli::cli_alert_info(
      "Cumulative variance at {.val {final_dims}} PCs: {.val {round(cumulative_variance[final_dims], 1)}}%"
    )
    cli::cli_alert_info(
      "Variance explained by PC{final_dims}: {.val {round(pct_variance[final_dims], 2)}}%"
    )

    # if (cumulative_variance[final_dims] < 80) {
    #     cli::cli_warn(
    #         "Cumulative variance < 80%. Consider increasing dimensions."
    #     )
    # }
  }

  final_dims
}
