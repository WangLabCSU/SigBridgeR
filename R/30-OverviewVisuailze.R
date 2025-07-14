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
FetchUMAP = function(
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
                    label = label,
                    order = order,
                    label.size = label_size,
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
        grid_size <- ifelse(
            plot_count < 4,
            plot_count,
            ceiling(sqrt(plot_count))
        )
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

#' @title Visualization of Cell Screening Fractions
#'
#' @description
#' Generates stacked bar plots showing the proportion of cells classified as
#' Positive/Negative/Neutral by single-cell screening algorithms (Scissor, scPAS,
#' scPP, or scAB) across different sample groups.
#'
#' @usage
#' ScreenFractionPlot(
#'   screened_seurat,
#'   group_by = "Source",
#'   screen_type = c("scissor", "scPAS", "scPP", "scAB"),
#'   show_null = FALSE,
#'   plot_color = NULL,
#'   show_plot = TRUE,
#'   plot_title = "Screen Fraction",
#'   stack_width = 0.85,
#'   x_text_angle = 45,
#'   axis_linewidth = 0.8,
#'   legend_position = "right"
#' )
#'
#' @param screened_seurat A Seurat object containing screening results in metadata.
#'        Must contain columns corresponding to `screen_type`.
#' @param group_by Metadata column name for grouping samples (default: "Source").
#' @param screen_type Screening algorithm used. Must match a metadata column
#'        (case-sensitive, e.g., "scissor" for Scissor results).
#' @param show_null Logical whether to show groups with zero cells (default: FALSE).
#' @param plot_color Custom color palette (named vector format):
#'        - Required names: "Positive", "Negative", "Neutral"
#'        - Default: c("Neutral"="#CECECE", "Positive"="#ff3333", "Negative"="#386c9b")
#' @param show_plot Logical to immediately display plot (default: TRUE).
#' @param plot_title Plot title (default: "Screen Fraction").
#' @param stack_width Bar width (default: 0.85).
#' @param x_text_angle X-axis label angle (default: 45).
#' @param axis_linewidth Axis line thickness (default: 0.8).
#' @param legend_position Legend position (default: "right").
#' @param x_lab X-axis label (default: NULL).
#' @param y_lab Y-axis label (default: "Status Fraction")
#'
#' @return A list containing:
#' \itemize{
#'   \item{stats: A data frame with screening statistics including:
#'     \itemize{
#'       \item Grouping variable counts
#'       \item Raw cell counts
#'       \item Percentage fractions
#'     }
#'   }
#'   \item{plot: A ggplot2 object of the stacked bar plot}
#' }
#'
#' @section Visualization Details:
#' - Bars are ordered by descending Positive fraction
#' - Y-axis shows percentage (0-100%)
#' - Zero-fraction groups are automatically hidden unless `show_null=TRUE`
#'
#' @examples
#' \dontrun{
#' # Basic usage with Scissor results
#' res <- ScreenFractionPlot(
#'   screened_seurat = scissor_result,
#'   group_by = "PatientID",
#'   screen_type = "scissor"
#' )
#'
#' # Customized plot
#' custom_plot <- ScreenFractionPlot(
#'   screened_seurat = scAB_result,
#'   group_by = "TissueType",
#'   screen_type = "scAB",
#'   plot_color = c("Positive"="red", "Negative"="blue", "Neutral"="grey80"),
#'   plot_title = "scAB Screening Results",
#'   x_text_angle = 90
#' )
#' }
#'
#' @export
#' @importFrom dplyr count group_by mutate ungroup pull filter arrange
#' @importFrom tidyr complete
#' @importFrom ggplot2 ggplot aes geom_col scale_y_continuous scale_fill_manual
#' @importFrom ggplot2 theme_classic labs theme element_text element_line
#' @importFrom scales percent_format
#' @importFrom glue glue
#'
#'
ScreenFractionPlot = function(
    screened_seurat,
    group_by = "Source",
    screen_type = c("scissor", "scPAS", "scPP", "scAB"),
    show_null = FALSE,
    plot_color = NULL,
    show_plot = TRUE,
    plot_title = "Screen Fraction",
    stack_width = 0.85,
    x_text_angle = 45,
    axis_linewidth = 0.8,
    legend_position = "right",
    x_lab = NULL,
    y_lab = "Status fraction"
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
                "scissor$|scPAS$|scPP$|scAB.*$",
                screened_seurat$scRNA_data@meta.data %>% names(),
                value = T
            ),
            .sep = "\n"
        ))
    }

    plot_color <- plot_color %||%
        c("Neutral" = "#CECECE", "Positive" = "#ff3333", "Negative" = "#386c9b")

    result <- list() # return

    stats_df <- screened_seurat@meta.data %>%
        dplyr::count(!!sym(group_by), !!sym(screen_type)) %>%
        tidyr::complete(
            !!sym(group_by),
            !!sym(screen_type),
            fill = list(n = 0)
        ) %>%
        dplyr::group_by(!!sym(group_by)) %>%
        dplyr::mutate(Total = sum(n)) %>% # total number of each group
        dplyr::ungroup() %>%
        dplyr::mutate(
            Fraction = ifelse(Total == 0, 0, n / Total)
        ) %>%
        {
            result$stats <<- .
            .
        } %>%
        dplyr::select(!!sym(group_by), !!sym(screen_type), Fraction)

    plot_order <- stats_df %>%
        dplyr::filter(!!sym(screen_type) == "Positive") %>%
        dplyr::arrange(desc(Fraction)) %>%
        dplyr::pull(!!sym(group_by))

    # for plot labels
    ms_type = screened_seurat@misc[[grep(
        screen_type,
        names(screened_seurat@misc),
        value = TRUE
    )]]

    # filter null records
    if (!show_null) {
        stats_df <- stats_df %>% dplyr::filter(Fraction > 0)
    }

    result$plot <- ggplot2::ggplot(
        stats_df,
        ggplot2::aes(
            x = factor(!!sym(group_by), levels = plot_order),
            y = `Fraction`,
            fill = !!sym(screen_type),
            title = plot_title
        )
    ) +
        ggplot2::geom_col(position = "stack", width = stack_width) +
        ggplot2::scale_y_continuous(
            labels = scales::percent_format(accuracy = 1),
            expand = c(0, 0),
            breaks = seq(0, 1, 0.1)
        ) +
        ggplot2::scale_fill_manual(
            values = plot_color
        ) +
        ggplot2::theme_classic(base_size = 14) +
        ggplot2::labs(
            x = x_lab,
            y = y_lab,
            fill = ms_type
        ) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(
                angle = x_text_angle,
                hjust = 1,
                vjust = 1
            ),
            axis.text = ggplot2::element_text(color = "black"),
            legend.position = legend_position,
            axis.line = ggplot2::element_line(linewidth = axis_linewidth)
        )

    if (show_plot) {
        print(result$plot)
    }

    return(result)
}
