# CellTypeAnnotationVisualize = function(seurat_obj, plot_color = NULL) {
#   if (
#     !requireNamespace(
#       c("Seurat", "dplyr", "cli", "SCpubr", "cowplot", "randomcoloR"),
#       quietly = TRUE
#     )
#   ) {
#     stop(
#       "Package `Seurat`, `dplyr`, `cli`, `cowplot`, `randomcoloR`  and `SCpubr` is required for this function"
#     )
#   }
#   if (!inherits(seurat_obj, "Seurat")) {
#     stop("Argument `seurat_obj` must be a Seurat object")
#   }
#   library(dplyr)

#   cli::cli_alert_info(c(
#     "[{TimeStamp()}] ",
#     crayon::green("Start Visualizing Cell Types of Seurat Object.")
#   ))
#   # cell types for annotation
#   cell_types <- unique(seurat_obj$cell_type)
#   cli::cli_alert_info("[{TimeStamp()}] Visualizing for:")
#   print(cell_types)
#   # colors of plots
#   plot_color = plot_color %||%
#     setNames(
#       randomcoloR::distinctColorPalette(k = length(cell_types), runTsne = TRUE),
#       cell_types
#     )
#   if (names(plot_color) |> is.null()) {
#     stop("`plot_color` must be a named vector.")
#   }

#   # basic umap
#   cli::cli_alert_info("[{TimeStamp()}] Visualizing for cell types")
#   cell_type_umap <- SCpubr::do_DimPlot(
#     sample = seurat_obj,
#     group.by = "cell_type",
#     label = TRUE,
#     legend.position = "bottom",
#     pt.size = 1.0,
#     label.size = 4,
#     label.box = TRUE,
#     repel = TRUE,
#     colors.use = plot_color,
#     plot.title = "Cell Type"
#   ) +
#     ggplot2::theme(
#       plot.title = ggplot2::element_text(
#         hjust = 0.5,
#         margin = ggplot2::margin(b = 15, t = 10)
#       ),
#       legend.text = ggplot2::element_text(size = 8),
#       legend.key.size = ggplot2::unit(0.3, "cm"),
#       plot.margin = ggplot2::unit(c(0.8, 0.8, 0.8, 0.8), "cm")
#     )
#   ifelse(
#     "consensus_proportion" %in% tolower(colnames(seurat_obj@meta.data)),
#     {
#       # Consensus Proportion Feature Plot
#       cli::cli_alert_info("[{TimeStamp()}] Visualizing for features")
#       consensus_proportion <- SCpubr::do_FeaturePlot(
#         sample = seurat_obj,
#         features = "consensus_proportion",
#         order = TRUE,
#         pt.size = 1.0,
#         enforce_symmetry = FALSE,
#         legend.title = "Consensus",
#         plot.title = "Consensus Proportion",
#         sequential.palette = "YlGnBu", # nature method standard
#         sequential.direction = 1, # from shallow to deep
#         min.cutoff = min(seurat_obj$consensus_proportion),
#         max.cutoff = max(seurat_obj$consensus_proportion),
#         na.value = "lightgrey"
#       ) +
#         ggplot2::theme(
#           plot.title = ggplot2::element_text(
#             hjust = 0.5,
#             margin = ggplot2::margin(b = 15, t = 10)
#           ),
#           plot.margin = ggplot2::unit(c(0.8, 0.8, 0.8, 0.8), "cm")
#         )

#       # Shannon Entropy Feature Plot
#       cli::cli_alert_info("[{TimeStamp()}] Visualizing for shannon features")
#       shannon_entropy <- SCpubr::do_FeaturePlot(
#         sample = seurat_obj,
#         features = "entropy",
#         order = TRUE,
#         pt.size = 1.0,
#         enforce_symmetry = FALSE,
#         legend.title = "Entropy",
#         plot.title = "Shannon Entropy",
#         sequential.palette = "OrRd", # nature method standard
#         sequential.direction = -1, # from deep to shallow
#         min.cutoff = min(seurat_obj$entropy),
#         max.cutoff = max(seurat_obj$entropy),
#         na.value = "lightgrey"
#       ) +
#         ggplot2::theme(
#           plot.title = ggplot2::element_text(
#             hjust = 0.5,
#             margin = ggplot2::margin(b = 15, t = 10)
#           ),
#           plot.margin = ggplot2::unit(c(0.8, 0.8, 0.8, 0.8), "cm")
#         )

#       # Combine Plots
#       combined_plot <- cowplot::plot_grid(
#         cell_type_umap,
#         consensus_proportion,
#         shannon_entropy,
#         ncol = 3,
#         rel_widths = c(1.2, 1.2, 1.2)
#       )
#     },
#     {
#       consensus_proportion = NULL
#       shannon_entropy = NULL
#       combined_plot = cell_type_umap
#     }
#   )

#   cli::cli_alert_info(
#     "[{TimeStamp()}] Visualization Finished"
#   )

#   return(list(
#     combined_plot = combined_plot,
#     cell_type_umap = cell_type_umap,
#     consensus_proportion = consensus_proportion,
#     shannon_entropy = shannon_entropy,
#     cell_types = cell_types
#   ))
# }

# #Visualization
# I recommend to save this plot to a file to see what it is like.
# score_heatmap = SingleR::plotScoreHeatmap(
#   prediction_sce,
#   order.by = "labels",
#   annotation_col = SummarizedExperiment::as.data.frame(SingleCellExperiment::colData(
#     sce_obj
#   )[, c("scissor", "Source", "Tissue"), drop = FALSE])
# )

# delta_distribution = SingleR::plotDeltaDistribution(
#   prediction_sce,
#   ncol = unique(prediction_sce$labels) |>
#     length() |>
#     sqrt() |>
#     ceiling() +
#     2
# )
