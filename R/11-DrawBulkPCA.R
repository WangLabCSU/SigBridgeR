#' Draw PCA Plot for Bulk Expression Data
#'
#' @description
#' Performs PCA on transposed bulk expression data and visualizes the first two
#' principal components as a scatter plot with marginal density plots. Optionally
#' overlays batch information as point shapes (up to 4 batches).
#'
#' The plot is composed via [patchwork] from:
#' - A top marginal density plot for PC1.
#' - A right marginal density plot for PC2.
#' - A central scatter plot with group-colored ellipses (via [ggforce::geom_mark_ellipse()]).
#'
#' @param bulk A numeric matrix or data frame of bulk expression data
#'   (genes x samples). Column names are used as sample identifiers.
#' @param group A vector of group labels, one per sample (column of `bulk`).
#'   Used for coloring points and ellipses.
#' @param batch An optional vector of batch labels, one per sample. When
#'   provided, points are additionally shaped by batch (up to 4 batches
#'   supported). Default: `NULL`.
#' @param show_plot Logical. If `TRUE` (default), the plot is printed
#'   immediately.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the combined `patchwork`/`ggplot` object.
#'
#' @examples
#' \dontrun{
#' # Small example matrix: 10 genes x 8 samples
#' set.seed(123)
#' bulk <- matrix(
#'   rnorm(100 * 80, mean = 10, sd = 2),
#'   nrow = 100, ncol = 80,
#'   dimnames = list(
#'     paste0("Gene", 1:100),
#'     paste0("Sample", 1:80)
#'   )
#' )
#'
#' group <- rep(c("Control", "Treatment"), each = 40)
#'
#' # Basic PCA plot without batch
#' DrawBulkPCA(bulk = bulk, group = group)
#'
#' # PCA plot with batch shapes
#' batch <- rep(c("Batch1", "Batch2"), times = 40)
#' DrawBulkPCA(bulk = bulk, group = group, batch = batch)
#' }
#'
#' @family input_preprocess
#' @export
DrawBulkPCA <- function(
  bulk,
  group,
  batch = NULL,
  show_plot = TRUE,
  ...
) {
  check_installed(c("ggplot2", "ggforce", "patchwork", "tibble"))

  pca <- stats::prcomp(t(bulk), scale. = TRUE)
  percent_var <- pca$sdev^2 / sum(pca$sdev^2)

  pca_df <- tibble::tibble(
    sample = colnames(bulk),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    group = group,
    batch = batch
  )

  PC1 <- PC2 <- NULL
  p_pca <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(x = `PC1`, y = `PC2`, color = `group`)
  ) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::labs(
      x = paste0("PC1 (", round(percent_var[1], 2), "%)"),
      y = paste0("PC2 (", round(percent_var[2], 2), "%)")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = c(0.95, 0.95), # 图例在右上角内部
      legend.justification = c(1, 1) # 对齐到右上角
    ) +
    ggforce::geom_mark_ellipse(
      ggplot2::aes(fill = group, group = group),
      alpha = 0.1,
      expand = ggplot2::unit(3, "mm"),
      show.legend = FALSE
    ) +
    ggplot2::geom_hline(yintercept = 0, color = "gray70", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, color = "gray70", linetype = "dashed")

  if ("batch" %chin% colnames(pca_df)) {
    n_batch <- length(unique(pca_df$batch))

    chk::chk_lt(n_batch, 5)

    p_pca <- p_pca +
      ggplot2::aes(shape = `batch`) +
      ggplot2::scale_shape_manual(
        values = c(
          16,
          17,
          18,
          19,
          20
        )[seq_len(n_batch)]
      )
  }

  pc1_range <- range(pca_df$PC1)
  pc2_range <- range(pca_df$PC2)

  # 创建密度图 - 描边颜色与填充颜色一致
  density_x <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(x = PC1, fill = group, color = group)
  ) +
    ggplot2::geom_density(alpha = 0.7, bw = "nrd", adjust = 2) +
    ggplot2::coord_cartesian(xlim = pc1_range) + # 与主图x轴范围一致
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")

  density_y <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(x = PC2, fill = group, color = group)
  ) +
    ggplot2::geom_density(alpha = 0.7, trim = FALSE, bw = "nrd", adjust = 1) +
    ggplot2::coord_cartesian(xlim = pc2_range) + # 与主图y轴范围一致
    ggplot2::coord_flip() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")

  combined_ellipse <- density_x +
    patchwork::plot_spacer() +
    p_pca +
    density_y +
    patchwork::plot_layout(ncol = 2, widths = c(4, 1), heights = c(1, 4)) +
    patchwork::plot_annotation(title = "Principal Component Analysis (PCA)")

  if (show_plot) {
    print(combined_ellipse)
  }

  invisible(combined_ellipse)
}
