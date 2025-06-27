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

# ---- Scissor(+/-/N)-sample/source stacked graph(in tumor) ----

SampleScissorFractionStackPlot = function(
  scissored,
  sample_colname = "Source",
  return_stats = FALSE,
  return_plot = TRUE,
  show_null = FALSE,
  plot_color = NULL,
  show_plot = TRUE
) {
  require(dplyr, quietly = TRUE)

  if (!inherits(scissored, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  if (!"scissor" %in% colnames(scissored@meta.data)) {
    stop("Missing 'scissor' metadata column")
  }
  if (!return_plot && !return_stats) {
    stop("This function is doing nothing now")
  }

  if (
    is.null(plot_color) ||
      length(plot_color) < 3 ||
      length(names(plot_color)) == 0
  ) {
    cli::cli_alert_warning(crayon::yellow(
      "WARNING: Color set not enough/wrong/not a named vector, now set to default colors."
    ))
    plot_color = stats::setNames(
      c("#CECECE", "#ff3333", "#386c9b"),
      c("N", "Pos", "Neg")
    )
  }
  
  stats_df <- scissored@meta.data %>%
    dplyr::count(!!sym(sample_colname), scissor) %>%
    tidyr::complete(
      !!sym(sample_colname),
      scissor = c("Neutral", "Positive", "Negative"),
      fill = list(n = 0)
    ) %>%
    dplyr::group_by(!!sym(sample_colname)) %>%
    dplyr::mutate(
      Total = sum(n),
      Status = factor(
        scissor,
        levels = c("Positive", "Negative", "Neutral"),
        labels = c("Pos", "Neg", "N")
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Fraction = ifelse(Total == 0, 0, n / Total),
      .keep = "unused"
    ) %>%
    dplyr::select(-scissor) %>%
    tidyr::pivot_wider(
      names_from = Status,
      values_from = Fraction,
      values_fill = 0
    ) %>%
    tidyr::pivot_longer(
      cols = c(Pos, Neg, N),
      names_to = "Scissor status",
      values_to = "Fraction"
    )

  # filter null records
  if (!show_null) {
    stats_df <- stats_df %>% dplyr::filter(Fraction > 0)
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
        values = stats::setNames(plot_color, c("N", "Pos", "Neg"))
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

# ---- 4. One_cell_type_fraction-tumor_in_Scissor(+)_fraction scatter plot ----

CellTypeMSScatterBoxplot = function(
  scissored,
  raw_data,
  cell_type_select = "T cell",
  ms_select = NULL,
  scissor_select = "Pos",
  plot_show = TRUE,
  plot_color = NULL,
  return_stats = FALSE,
  return_plot = TRUE
) {
  if (!inherits(scissored, "Seurat") || !inherits(raw_data, "Seurat")) {
    stop("Both inputs  `scissored` and `raw_data` must be Seurat objects")
  }

  # factor-scissor_select transfer
  scissor_type <- dplyr::case_when(
    tolower(scissor_select) %in% c("pos", "positive", "p") ~ "Positive",
    tolower(scissor_select) %in% c("neg", "negative") ~ "Negative",
    TRUE ~ "Neutral"
  )

  meta_raw <- raw_data@meta.data
  meta_scissor <- scissored@meta.data

  x_data <- meta_raw %>%
    dplyr::group_by(Source) %>%
    dplyr::summarise(
      Tumor_count = sum(cnv_status == "tumor"),
      Immune_count = sum(
        Tissue
      ),
      Celltype_count = sum(Tissue == cell_type_select),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      meta_scissor %>%
        dplyr::group_by(Source) %>%
        dplyr::summarise(
          Scissored_tumor = sum(
            scissor == scissor_type # these cells are tumor cells
          ),
        ),
      by = "Sample"
    ) %>%
    dplyr::filter(dplyr::if_all(dplyr::everything(), ~ !is.na(.) & . != 0)) %>%
    dplyr::mutate(
      Tumor_fraction = Scissored_tumor / Tumor_count,
      Immune_fraction = Celltype_count / Immune_count,
      Tumor_fraction = ifelse(is.na(Tumor_fraction), 0, Tumor_fraction)
    )

  # linear regression
  model <- lm(Immune_fraction ~ Tumor_fraction, data = x_data)
  model_summary <- summary(model)

  # plot color
  if (plot_show || return_plot) {
    plot_color <- plot_color %||%
      c("#ff3333", "#660000", "black", '#E2C8AF', '#B6DDEB')
    # scatter plot
    gg_scatter <- ggplot2::ggplot(
      x_data,
      ggplot2::aes(x = Tumor_fraction, y = Immune_fraction)
    ) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method = "lm", se = FALSE, color = plot_color[1]) +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = Sample),
        box.padding = 0.35,
        point.padding = 0.5,
        size = 3,
        color = plot_color[3]
      ) +
      ggplot2::annotate(
        "text",
        x = min(x_data$Tumor_fraction, na.rm = TRUE) +
          0.1 * diff(range(x_data$Tumor_fraction, na.rm = TRUE)),
        y = max(x_data$Immune_fraction, na.rm = TRUE) -
          0.1 * diff(range(x_data$Immune_fraction, na.rm = TRUE)),
        label = sprintf(
          "y = %.4f %s %.4fx\nr = %s%.4f\np = %.4f",
          coef(model)[1],
          ifelse(coef(model)[2] >= 0, "+", "-"),
          abs(coef(model)[2]),
          ifelse(coef(model)[2] >= 0, "", "-"),
          abs(sqrt(model_summary$r.squared)),
          model_summary$coefficients[2, 4]
        ),
        color = plot_color[1],
        hjust = 0,
        size = 4
      ) +
      ggplot2::scale_x_continuous(labels = scales::percent) +
      ggplot2::scale_y_continuous(labels = scales::percent) +
      cowplot::theme_cowplot(12) +
      ggplot2::labs(
        x = "Scissored tumor cell fraction",
        y = glue::glue("{cell_type_select} fraction")
      )

    # boxplot
    x_data <- x_data %>%
      dplyr::mutate(
        Group = ifelse(
          Tumor_fraction <= median(Tumor_fraction),
          "BotHalf",
          "TopHalf"
        )
      )
    gg_box <- ggplot2::ggplot(
      x_data,
      ggplot2::aes(x = Group, y = Immune_fraction, fill = Group)
    ) +
      ggplot2::geom_boxplot() +
      ggplot2::stat_summary(
        fun = "max",
        geom = "errorbar",
        aes(ymax = ..y.., ymin = ..y..),
        width = 0.5,
        color = "black"
      ) +
      ggplot2::stat_summary(
        fun = "min",
        geom = "errorbar",
        aes(ymax = ..y.., ymin = ..y..),
        width = 0.5,
        color = "black"
      ) +
      ggplot2::geom_jitter(width = 0.2) +
      ggplot2::scale_fill_manual(values = plot_color[4:5]) +
      cowplot::theme_cowplot(12) +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = Sample),
        box.padding = 0.35,
        point.padding = 0.5,
        size = 3,
        color = plot_color[3]
      ) +
      ggplot2::labs(y = glue::glue("{cell_type_select} fraction"))

    if (plot_show) {
      print(patchwork::wrap_plots(gg_box, gg_scatter))
    }
  }

  if (return_stats) {
    wilcox_test <- wilcox.test(Immune_fraction ~ Group, data = x_data)

    stats_all = list(
      Cell_type = cell_type_select,
      R_squared = model_summary$r.squared,
      P_value_lm = model_summary$coefficients[2, 4],
      P_value_wilcox = wilcox_test$p.value,
      ms_select = ms_select
    )
  }

  result = list()
  if (return_plot) {
    result$scatter_plot = gg_scatter
    result$box_plot = gg_box
  }
  if (return_stats) {
    result$stats =
      stats_all
  }
  return(result)
}

DoAllCelltypeMS <- function(
  scissored,
  raw_data,
  scissor_select = "p",
  ms_select,
  plot_show = FALSE,
  plot_color = NULL,
  plot_save = TRUE,
  return_stats = TRUE,
  return_plot = TRUE,
  plot_save_dir = "./plot"
) {
  require(cli, quietly = T)
  cli::cli_alert_info(
    "Analysis is being conducted on {crayon::green('all types of immune cells')}."
  )
  cell_types <- levels(raw_data@meta.data$Celltype)
  cell_types = cell_types[
    !cell_types %in%
      c(
        "Endothelial",
        "Epithelial",
        "Fibroblast"
      )
  ]
  if (is.null(cell_types)) {
    stop("No available immune cell types ")
  }

  stats_data = do.call(
    rbind,
    lapply(cell_types, function(celltype) {
      res = tryCatch(
        {
          record <- CellTypeMSScatterBoxplot(
            scissored = scissored,
            raw_data = raw_data,
            cell_type_select = celltype,
            scissor_select = scissor_select,
            plot_show = plot_show,
            plot_color = plot_color,
            return_stats = return_stats,
            return_plot = return_plot,
            ms_select = ms_select
          )
          if (plot_save) {
            ggsave(
              file.path(
                plot_save_dir,
                glue::glue("figure3_{celltype}_tumor_fraction_scatter.pdf")
              ),
              plot = record$scatter_plot,
              width = 10,
              height = 8,
              dpi = 300
            )
            ggsave(
              file.path(
                plot_save_dir,
                glue::glue("figure3_{celltype}_tumor_fraction_box.pdf")
              ),
              plot = record$box_plot,
              width = 10,
              height = 8,
              dpi = 300
            )
          }
          cli::cli_alert_success(glue::glue(
            "{format(Sys.time(), '%Y/%m/%d %H:%M:%S')} {crayon::green('Success:')} {celltype}"
          ))

          return(invisible(record$stats))
        },
        error = function(e) {
          cli::cli_alert_danger(glue::glue("Error processing: {e$message}"))
          return(data.frame("error"))
        }
      )
      if (is.null(res)) data.frame() else res
    })
  )
  return(stats_data)
}

# from fig3
DoAllCelltypeMS.Heatmap = function(
  all_stats_data,
  title = NULL,
  param_type = c("P_value_lm", "P_value_wilcox", "R_squared"),
  color_spectrum = NULL,
  plot_show = FALSE,
  rowname_fontsize = 10,
  colname_fontsize = 10,
  group_label_fontsize = 8,
  heatmap_width = 12,
  heatmap_height = 8
) {
  RequireNamespace <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(glue::glue("Package {crayon::red(pkg)} required but not installed"))
    }
  }
  RequireNamespace("ComplexHeatmap")
  RequireNamespace("data.table")
  RequireNamespace("circlize")
  library(dplyr)
  if (is.null(color_spectrum)) {
    color_spectrum = circlize::colorRamp2(
      c(0, 0.5, 1),
      c("#117011", "white", "#FF0000")
    )
  }
  if (!all(param_type %in% c("P_value_lm", "P_value_wilcox", "R_squared"))) {
    stop(glue::glue(
      "{crayon::red('param_type')} set wrong.",
      "valid 'param_type': P_value_lm, P_value_wilcox, R_squared",
      .sep = "\n"
    ))
  }

  all_ms_type <- names(all_stats_data)
  all_cell_types <- unique(unlist(lapply(
    all_stats_data,
    function(x) {
      if ("Cell_type" %in% colnames(x)) {
        as.character(x[, "Cell_type"])
      } else {
        cli::cli_alert_warning(
          "Column 'Cell_type' missing in sample: {crayon::yellow(names(x)[1])} "
        )
        character(0)
      }
    }
  )))

  cli::cli_alert_info(glue::glue(
    "Drawing heatmap for {crayon::bold(length(all_ms_type))} samples",
    "{crayon::bold(length(all_cell_types))} cell types",
    .sep = ", "
  ))

  mat <- data.table::rbindlist(
    lapply(all_ms_type, function(ms) {
      stats <- all_stats_data[[ms]]
      data.table::data.table(
        Cell_type = unlist(stats[, 1]),
        Mutational_signature = ms,
        R_squared = unlist(stats[, 2]),
        P_value_lm = unlist(stats[, 3]),
        P_value_wilcox = unlist(stats[, 4])
      )
    }),
    fill = TRUE
  ) %>%
    data.table::melt(
      id.vars = c("Cell_type", "Mutational_signature"),
      measure.vars = param_type,
      variable.name = "Parameter",
      value.name = "Value"
    ) %>%
    dplyr::mutate(Value = suppressWarnings(as.numeric(Value))) %>%
    data.table::dcast(
      Cell_type ~ Mutational_signature + Parameter,
      value.var = "Value"
    ) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("Cell_type") %>%
    as.matrix() %>%
    `mode<-`("numeric")

  GeneratGroupColors <- function(n_groups) {
    if (n_groups == 1) {
      return("#F16913")
    }
    warm_pal <- grDevices::colorRampPalette(c(
      "#DFD3E5",
      "#CCB8D7",
      "#C6A2D9",
      "#AF8EC0",
      "#997AA9",
      "#8A6E99"
    ))
    cool_pal <- grDevices::colorRampPalette(c(
      "#EFF3FF",
      "#C6DBEF",
      "#9ECAE1",
      "#6BAED6",
      "#4292C6",
      "#2171B5"
    ))
    n_warm <- ceiling(n_groups / 2)
    n_cool <- floor(n_groups / 2)
    return(c(
      warm_pal(n_warm)[seq_len(n_warm)],
      cool_pal(n_cool)[seq_len(n_cool)]
    )[order(c(seq(1, n_groups, 2), seq(2, n_groups, 2)))])
  }

  ht <- ComplexHeatmap::Heatmap(
    matrix = mat,
    col = color_spectrum,
    column_labels = gsub(".*_(lm|wilcox|R_squared)$", "\\1", colnames(mat)),
    column_split = rep(all_ms_type, each = length(param_type)),
    column_names_gp = grid::gpar(fontsize = colname_fontsize),
    row_title_gp = grid::gpar(fontsize = rowname_fontsize),
    row_order = rownames(mat),
    na_col = "grey90",
    rect_gp = grid::gpar(col = "white", lwd = 1),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      Sample = ComplexHeatmap::anno_block(
        gp = grid::gpar(
          col = "white",
          fill = GeneratGroupColors(length(all_ms_type))
        ),
        labels = all_ms_type,
        labels_gp = grid::gpar(
          fontsize = group_label_fontsize,
          col = "black",
          fontface = "bold"
        )
      )
    ),
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(title = "Value"),
    column_title = title, # plot title
    row_title = "Cell Types",
    cluster_columns = FALSE,
    heatmap_width = grid::unit(heatmap_width, "cm"),
    heatmap_height = grid::unit(heatmap_height, "cm")
  )

  cli::cli_alert_success(glue::glue(
    "{crayon::green('Success:')} heatmap visualizatioin."
  ))
  if (plot_show) {
    ComplexHeatmap::draw(ht)
  }
  return(ht)
}
