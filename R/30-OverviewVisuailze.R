# ----  Visualization ----

#' @title Generate UMAP Plots for Seurat Object
#'
#' @description
#' Creates customizable UMAP visualizations for single-cell data, supporting both
#' cluster visualization and gene expression plotting. Handles multiple plotting
#' parameters through flexible argument passing.
#'
#' @usage
#' FetchUAMP(
#'   seurat_obj,
#'   group_by = c("celltype"),
#'   feature = NULL,
#'   plot_name = NULL,
#'   plot_show = FALSE,
#'   point_size = 0.6,
#'   plot_color = NULL,
#'   label_size = 10,
#'   umap_label = TRUE,
#'   order = c(2, 1),
#'   feature_max_cutoff = 2,
#'   feature_min_cutoff = -2,
#'   feature_plot_raster = FALSE,
#'   ...
#' )
#'
#' @param seurat_obj A Seurat object containing UMAP reduction
#' @param group_by Metadata column(s) for cluster visualization (vector).
#'        Set to NULL to disable.
#' @param feature Feature(s) to plot expression (vector). Set to NULL to disable.
#' @param plot_name Base title for cluster plots
#' @param plot_show Logical to display combined plot immediately (default: FALSE)
#' @param point_size Point size for UMAP (default: 0.6)
#' @param plot_color Color palette for clusters (default: Seurat palette)
#' @param label_size Label font size (default: 10)
#' @param label Logical to show cluster labels (default: TRUE)
#' @param order Plotting order for points (default: c(2,1) - Significant first)
#' @param feature_max_cutoff Maximum expression cutoff (default: NA)
#' @param feature_min_cutoff Minimum expression cutoff (default: NA)
#' @param feature_plot_raster Logical to rasterize feature plots (default: NULL)
#' @param ... Additional arguments passed to either:
#'        - `Seurat::DimPlot()` for cluster plots
#'        - `Seurat::FeaturePlot()` for expression plots
#'
#' @return A list of ggplot objects containing:
#' \itemize{
#'   \item For `group_by`: UMAP cluster plots
#'   \item For `feature`: Feature expression plots
#' }
#'
#' @section Plotting Logic:
#' 1. For each column in `group_by`, creates a separate DimPlot
#' 2. For each entry in `feature`, creates a separate FeaturePlot
#' 3. Combines all plots when `plot_show=TRUE` using patchwork
#'
#' @examples
#' \dontrun{
#' # Cluster visualization only
#' plots <- FetchUAMP(
#'   seurat_obj = pbmc,
#'   group_by = c("celltype", "cluster"),
#'   plot_name = "PBMC"
#' )
#'
#' # Feature expression with custom parameters
#' plots <- FetchUAMP(
#'   seurat_obj = pbmc,
#'   feature = c("CD3D", "CD8A"),
#'   feature_max_cutoff = 3,
#'   feature_min_cutoff = 0,
#' )
#'
#' }
#'
#' @export
#' @importFrom Seurat DimPlot FeaturePlot
#' @importFrom ggplot2 ggtitle
#' @importFrom patchwork wrap_plots
#' @importFrom glue glue
#'
FetchUAMP = function(
  seurat_obj,
  group_by = NULL,
  feature = NULL,
  plot_name = NULL,
  plot_show = FALSE,
  point_size = 0.6,
  plot_color = NULL,
  label_size = 10,
  label = TRUE,
  order = c(2, 1),
  feature_max_cutoff = NA,
  feature_min_cutoff = NA,
  feature_plot_raster = NULL,
  ...
) {
  extra_params <- list(...)

  dimplot_args <- names(formals(Seurat::DimPlot))
  featureplot_args <- names(formals(Seurat::FeaturePlot))

  dimplot_dot <- extra_params[names(extra_params) %in% dimplot_args]
  featureplot_dot <- extra_params[names(extra_params) %in% featureplot_args]

  umap_list <- list()

  if (!is.null(group_by)) {
    umap_list <- lapply(X = group_by, FUN = function(group_name) {
      dim_args <- c(
        list(
          object = seurat_obj,
          reduction = "umap",
          label = label,
          group.by = group_name,
          cols = plot_color,
          pt.size = point_size,
          order = order,
          label.size = label_size
        ),
        dimplot_dot
      )

      umap_plot <- do.call(Seurat::DimPlot, dim_args) +
        ggplot2::ggtitle(glue::glue("{plot_name} UMAP"))
      return(umap_plot)
    })
  }

  if (!is.null(feature)) {
    feature_plots <- lapply(feature, function(feat) {
      feat_args <- c(
        list(
          object = seurat_obj,
          features = feat,
          cols = plot_color,
          raster = feature_plot_raster,
          label = label,
          order = order,
          label.size = label_size,
          max.cutoff = feature_max_cutoff,
          min.cutoff = feature_min_cutoff
        ),
        featureplot_dot
      )

      feature_plot <- do.call(Seurat::FeaturePlot, feat_args) +
        ggplot2::ggtitle(glue::glue("{feat} Expression"))

      return(feature_plot)
    })
    umap_list <- c(umap_list, feature_plots)
  }

  plot_count <- length(umap_list)
  if (plot_count == 0) {
    stop("No plots have been generated")
  }
  names(umap_list) <- c(group_by, feature)

  if (plot_show) {
    grid_size <- ifelse(plot_count < 4, plot_count, ceiling(sqrt(plot_count)))
    combined_plot <- patchwork::wrap_plots(
      umap_list,
      ncol = grid_size,
      nrow = ceiling(plot_count / grid_size)
    )
    print(combined_plot)
  }
  return(umap_list)
}

# ---- Screened cell fraction(+/-/N)-sample/source stacked graph ----

ScreenFractionPlot = function(
  screened_seurat,
  group_by = "Source",
  screen_type = c("scissor", "scPAS", "scPP", "scAB"),
  return_stats = TRUE,
  return_plot = TRUE,
  show_null = FALSE,
  plot_color = NULL,
  show_plot = TRUE
) {
  library(dplyr)

  if (!inherits(screened_seurat, "Seurat")) {
    stop("`screened_seurat` must be a Seurat object")
  }

  if (length(screen_type) != 1) {
    stop(glue::glue(
      "Please refer one screen algorithm type.",
      "Available screen types: ",
      grep(
        "scissor$|scPAS$|[Ss]cPP$|scAB.*$",
        screened_seurat$scRNA_data@meta.data %>% names(),
        value = T
      ),
      .sep = "\n"
    ))
  }

  plot_color <- plot_color %||%
    stats::setNames(
      c("Neutral", "Positive", "Negative"),
      c("#CECECE", "#ff3333", "#386c9b")
    )

  stats_df <- screened_seurat@meta.data %>%
    dplyr::count(!!sym(group_by), screen_type) %>%
    tidyr::complete(
      !!sym(group_by),
      screen_type = c("Neutral", "Positive", "Negative"),
      fill = list(n = 0)
    ) %>%
    dplyr::group_by(!!sym(group_by)) |>
    dplyr::mutate(
      Total = sum(n),
      Status = factor(
        screen_type,
        levels = c("Positive", "Negative", "Neutral"),
        labels = c("Positive", "Negative", "Neutral")
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      Fraction = ifelse(Total == 0, 0, n / Total),
      .keep = "unused"
    ) |>
    dplyr::select(-screen_type) |>
    tidyr::pivot_wider(
      names_from = Status,
      values_from = Fraction,
      values_fill = 0
    ) |>
    tidyr::pivot_longer(
      cols = c(Positive, Negative, Neutral),
      names_to = glue::glue("{screen_type} status"),
      values_to = "Fraction"
    )

  # filter null records
  if (!show_null) {
    stats_df <- stats_df |> dplyr::filter(Fraction > 0)
  }
  if (show_plot || return_plot) {
    gg <- ggplot2::ggplot(
      stats_df,
      ggplot2::aes(
        x = !!sym(group_by),
        y = `Fraction`,
        fill = `Scissor status`
      )
    ) +
      ggplot2::geom_col(position = "stack", width = 0.85) +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(accuracy = 1),
        expand = c(0, 0),
        breaks = seq(0, 1, 0.1)
      ) +
      ggplot2::scale_fill_manual(
        values = plot_color
      ) +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::labs(x = NULL, y = "Scissor status fraction") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = ggplot2::element_text(color = "black"),
        legend.position = "right",
        axis.line = ggplot2::element_line(linewidth = 0.8)
      )
    if (show_plot) {
      print(gg)
    }
  }

  result <- list()
  if (return_stats) {
    result$stats <- as.data.frame(stats_df)
  }
  if (return_plot) {
    result$plot <- gg
  }
  if (length(result) > 0) return(result)
}
