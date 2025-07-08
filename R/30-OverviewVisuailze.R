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
#' @param umap_label Logical to show cluster labels (default: TRUE)
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
#'   pt.size = 1.2
#' )
#'
#' # Combined cluster and feature plots
#' plots <- FetchUAMP(
#'   seurat_obj = pbmc,
#'   group_by = "celltype",
#'   feature = "CD79A",
#'   plot_show = TRUE
#' )
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
  umap_label = TRUE,
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
          label = umap_label,
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

  if (!is.null(feature)) {<<<HERE
    feature_plots <- lapply(feature, function(feat) {
      feat_args <- c(
        list(
          object = seurat_obj,
          features = feat,
          raster = feature_plot_raster,
          label = umap_label,
          order = order,
          label.size = label_size,
          max.cutoff = feature_max_cutoff,
          min.cutoff = feature_min_cutoff
        ),
        featureplot_dot
      )

      do.call(Seurat::FeaturePlot, feat_args) +
        ggplot2::ggtitle(glue::glue("{feat} Expression"))
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

# ---- The proportion of cell types in different samples/sources ----

CalculateCellTypeFraction = function(
  seuratobject,
  cell_type_col = "Celltype",
  select_cell_type = "All",
  return_stats = FALSE,
  return_plot = TRUE,
  plot_color = NULL,
  show_plot = TRUE,
  title = NULL
) {
  require(dplyr, quietly = TRUE)

  if (!inherits(seuratobject, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  cell_types <- seuratobject[[cell_type_col]] %>%
    unique() %>%
    unlist() %>%
    na.omit()

  # default: select all
  if (any(tolower(select_cell_type) %in% c("all", "", NULL))) {
    select_cell_type <- cell_types
  }
  # check if celltype inputed is wrong
  invalid_types <- setdiff(select_cell_type, cell_types)
  if (length(invalid_types) > 0) {
    cli::cli_alert_warning(glue::glue(
      "Available cell types:\n",
      glue::glue(cell_types, .sep = "\n")
    ))
    stop("Invalid cell types: ", glue::glue(invalid_types, collapse = ", "))
  }

  stats_df <- seuratobject@meta.data %>%
    dplyr::filter(!!rlang::sym(cell_type_col) %in% select_cell_type) %>%
    dplyr::count(Source, !!rlang::sym(cell_type_col), name = "Type_count") %>%
    dplyr::group_by(Source) %>%
    dplyr::mutate(
      Total_count = sum(Type_count),
      Fraction = Type_count / Total_count
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(Source, dplyr::desc(Fraction))

  if (return_plot || show_plot) {
    plot_color <- plot_color %||%
      randomcoloR::distinctColorPalette(
        length(unique(seuratobject@meta.data[[
          cell_type_col
        ]])),
        runTsne = TRUE
      )

    gg <- ggplot2::ggplot(
      stats_df,
      ggplot2::aes(
        x = `Source`,
        y = `Fraction`,
        fill = !!rlang::sym(cell_type_col)
      )
    ) +
      ggplot2::geom_col(position = ggplot2::position_stack(reverse = TRUE)) +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(),
        expand = c(0, 0)
      ) +
      ggplot2::scale_fill_manual(values = plot_color) +
      ggplot2::labs(
        x = NULL,
        y = "Cell Fraction",
        fill = cell_type_col,
        title = title
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.x = ggplot2::element_blank(),
        legend.position = "right",
        plot.margin = ggplot2::margin(1, 1, 1, 2, "cm"),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 15)
      )

    if (show_plot) {
      show(gg)
    }
  }

  return_list <- list()
  if (return_stats) {
    return_list$stats <- as.data.frame(stats_df)
  }
  if (return_plot) {
    return_list$plot <- gg
  }
  return(return_list)
}

# ---- Scissor(+/-/N)-sample/source stacked graph(in tumor) ----

SampleScreenFractionStackPlot = function(
  screened_seurat,
  sample_colname = "Source",
  screen_type = c("scissor", "scPAS", "scPP", "scAB"),
  return_stats = FALSE,
  return_plot = TRUE,
  show_null = FALSE,
  plot_color = NULL,
  show_plot = TRUE
) {
  require(dplyr, quietly = TRUE)

  if (!inherits(screened_seurat, "Seurat")) {
    stop("`screened_seurat` must be a Seurat object")
  }

  if (length(screen_type) != 1) {
    stop(glue::glue(
      "Please refer one screen algorithm type",
      "Available screen types: ",
      c('scissor', 'scPAS', 'scPP', 'scAB')[
        c('scissor', 'scPAS', 'scPP', 'scAB') %in%
          colnames(screened_seurat@meta.data)
      ],
      .sep = "\n"
    ))
  }

  if (!return_plot & !return_stats) {
    stop(
      "This function is doing nothing now, please set `return_stats` or `return_plot` to `TRUE`"
    )
  }

  plot_color <- plot_color %||%
    stats::setNames(c("N", "Pos", "Neg"), c("#CECECE", "#ff3333", "#386c9b"))

  stats_df <- screened_seurat@meta.data %>%
    dplyr::count(!!sym(sample_colname), screen_type) %>%
    tidyr::complete(
      !!sym(sample_colname),
      screen_type = c("Neutral", "Positive", "Negative"),
      fill = list(n = 0)
    ) %>%
    dplyr::group_by(!!sym(sample_colname)) |>
    dplyr::mutate(
      Total = sum(n),
      Status = factor(
        screen_type,
        levels = c("Positive", "Negative", "Neutral"),
        labels = c("Pos", "Neg", "N")
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
      cols = c(Pos, Neg, N),
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
        x = !!sym(sample_colname),
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
