##' @importFrom ggplot2 autoplot
##' @export
ggplot2::autoplot

#' autoplot  method for Seurat UMAP
#'
#' @description
#' A ggplot2 autoplot method for class ScreenMethodResult
#'
#'
#' @param object Seurat object
#' @param group.by Metadata column used for grouping/coloring
#' @param label Whether to label clusters
#' @param label.size Label font size
#' @param pt.size Point size
#' @param cols Custom colors
#' @param title Plot title
#' @param save_path If not NULL, save plot to this path
#' @param width Figure width
#' @param height Figure height
#' @param ... Other arguments passed to [Seurat::DimPlot]
#'
#' @return ggplot object
#' @export
`autoplot.SigBridgeR::ScreenMethodResult` <- function(
  object,
  group.by = NULL,
  label = TRUE,
  label.size = 4,
  pt.size = 0.5,
  cols = NULL,
  title = NULL,
  save_path = NULL,
  width = 6,
  height = 5,
  ...
) {
  check_installed(c("tidydr", "ggplot2"))

  cluster_color <- cols %||% palette_SigBridgeR()

  p <- Seurat::DimPlot(
    object = object@scRNA_data,
    reduction = "umap",
    label = label,
    label.size = label.size,
    group.by = group.by,
    cols = cluster_color,
    pt.size = pt.size,
    ...
  ) +
    Seurat::NoAxes() +
    ggplot2::labs(x = "UMAP_1", y = "UMAP_2") +
    tidydr::theme_dr() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::ggtitle(title)

  if (!is.null(save_path)) {
    ggplot2::ggsave(
      plot = p,
      filename = save_path,
      width = width,
      height = height,
      dpi = 400
    )
  }

  p
}

#' @title A Palette
#' @export
palette_SigBridgeR <- \() {
  c(
    "#B8D9A0",
    "#E8E0DB",
    "#A8D9D8",
    "#F5C0A8",
    "#D5C9E8",
    "#9BBBD9",
    "#B8A8E5",
    "#D5E5C0",
    "#E8D8E8",
    "#88D5B0",
    "#E5B8E5",
    "#88E0D5",
    "#A8D5A8",
    "#E5A8D0",
    "#F5C8E5",
    "#F0E890",
    "#F5E8B8",
    "#B8BCA0",
    "#C8E5D0",
    "#B888D5",
    "#C5B8E8",
    "#E5C090",
    "#D5A0A8",
    "#88A8D5",
    "#F5D890",
    "#D8E5B0",
    "#B8D5D0",
    "#E0E8D8",
    "#E0A8E5",
    "#A888D5",
    "#C5B0B8",
    "#E5A8E8",
    "#E8A8D5",
    "#D0C5B8",
    "#F5A098",
    "#B0E5D0",
    "#C0D8E5",
    "#90C8D5",
    "#E8A0B0",
    "#F5E0B8"
  )
}
