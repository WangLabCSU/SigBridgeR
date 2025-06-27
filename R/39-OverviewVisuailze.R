# ---- Tumor UMAP Visualization ----

FetchUMAP = function(
  SeuratObject,
  plot_name = NULL,
  plot_show = FALSE,
  group_by = NULL,
  point_size = 0.6,
  plot_color = NULL,
  label_size = 10,
  order = c(2, 1)
) {
  umap_plot = Seurat::DimPlot(
    object = SeuratObject,
    reduction = "umap",
    label = TRUE,
    group.by = group_by,
    cols = plot_color,
    pt.size = point_size,
    order = order,
    label.size = label_size
  ) +
    ggplot2::ggtitle(glue::glue("{plot_name} UMAP"))

  if (plot_show) {
    print(umap_plot)
  }
  return(umap_plot)
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
  if (any(tolower(select_cell_type) %in% c("all", ""))) {
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
