# ---- Visualization ----

# ---- Screened cell fraction(+/-/N)-sample/source stacked graph ----

#' @title Visualization of Cell Screening Fractions
#'
#' @description
#' Generates stacked bar plots showing the proportion of cells classified as
#' Positive/Negative/Neutral by single-cell screening algorithms (Scissor, scPAS,
#' scPP, or scAB) across different sample groups. Supports both single and multiple
#' screen types visualization.
#'
#' @usage
#' ScreenFractionPlot(
#'   screened_seurat,
#'   group_by = "Source", # x-axis label
#'   screen_type = NULL, # default: search all available screens
#'   show_null = FALSE,
#'   plot_color = NULL, # stack bar color of each status
#'   show_plot = FALSE,
#'   plot_title = "Screen Fraction",
#'   stack_width = 0.85,
#'   x_text_angle = 45,
#'   axis_linewidth = 0.8,
#'   x_lab = NULL,
#'   y_lab = "Fraction of Status",
#'   ncol = 2, # number of columns for facet wrap
#'   nrow = NULL, # number of rows for facet wrap
#'   verbose = SigBridgeRUtils::getFuncOption("verbose"),
#'   ... # ggplot2 theme arguments
#' )
#'
#' @param screened_seurat A Seurat object containing screening results in metadata.
#'        Must contain columns corresponding to `screen_type`.
#' @param group_by Metadata column name for grouping samples (default: "Source").
#' @param screen_type Screening algorithm(s) used. Can be a single value or vector.
#'        Must match metadata column(s) (case-sensitive, e.g., "scissor" for Scissor results).
#' @param show_null Logical whether to show groups with zero cells (default: FALSE).
#' @param plot_color Custom color palette (named vector format):
#'        - Required names: "Positive", "Negative", "Neutral"
#'        - Default: c("Neutral"="#CECECE", "Other"="#CECECE", Positive"="#ff3333", "Negative"="#386c9b")
#' @param show_plot Logical to immediately display plot (default: TRUE).
#' @param plot_title Plot title (default: "Screen Fraction"). When multiple screen types,
#'        can be a vector of titles or single title (will append screen type prefix).
#' @param stack_width Bar width (default: 0.85).
#' @param x_text_angle X-axis label angle (default: 45).
#' @param axis_linewidth Axis line thickness (default: 0.8).
#' @param x_lab X-axis label (default: NULL).
#' @param y_lab Y-axis label (default: "Fraction of Status").
#' @param ncol Number of columns for facet wrap when multiple screen types (default: 2).
#' @param nrow Number of rows for facet wrap when multiple screen types (default: NULL).
#' @param verbose Logical, whether to print a message
#' @param ... Other arguments passed to ggplot2::theme()
#'
#' @return A list containing:
#' \itemize{
#'   \item{stats: A data frame (single screen) or list of data frames (multiple screens)
#'     with screening statistics including:
#'     \itemize{
#'       \item Grouping variable counts
#'       \item Raw cell counts
#'       \item Percentage fractions
#'     }
#'   }
#'   \item{plot: A ggplot2 object (single screen) or list of ggplot2 objects (multiple screens)}
#'   \item{combined_plot: A combined plot using patchwork (only for multiple screens)}
#' }
#'
#' @section Visualization Details:
#' - Bars are ordered by descending Positive fraction
#' - Y-axis shows percentage (0-100%)
#' - Zero-fraction groups are automatically hidden unless `show_null=TRUE`
#' - For multiple screen types, plots can be combined using patchwork
#'
#' @examples
#' \dontrun{
#' # Single screen type usage
#' res <- ScreenFractionPlot(
#'   screened_seurat = scissor_result,
#'   group_by = "PatientID",
#'   screen_type = "scissor"
#' )
#'
#' # Multiple screen types at once
#' multi_res <- ScreenFractionPlot(
#'   screened_seurat = multi_screened_result,
#'   group_by = "TissueType",
#'   screen_type = c("scissor", "scPAS", "scPP", "scAB", "DEGAS"),
#'   ncol = 2
#' )
#' }
#'
#' @export
#' @family visualization_function
#'
ScreenFractionPlot <- function(
  screened_seurat,
  group_by = "Source",
  screen_type = NULL,
  show_null = FALSE,
  plot_color = NULL,
  show_plot = FALSE,
  plot_title = "Screen Fraction",
  stack_width = 0.85,
  x_text_angle = 45L,
  axis_linewidth = 0.8,
  x_lab = NULL,
  y_lab = "Fraction of Status",
  ncol = 2L,
  nrow = NULL,
  verbose = SigBridgeRUtils::getFuncOption("verbose"),
  ...
) {
  chk::chk_is(screened_seurat, "Seurat")
  chk::chk_character(group_by)
  chk::chk_length(group_by, 1)
  chk::chk_flag(show_null)
  if (!is.null(plot_color)) {
    chk::chk_vector(plot_color)
    chk::chk_named(plot_color)
  }
  check_installed(c("ggplot2", "patchwork", "dplyr"))

  dots <- rlang::list2(...)
  theme_args <- SigBridgeRUtils::FilterArgs4Func(dots, ggplot2::theme)

  meta_data <- screened_seurat[[]]
  all_screen_types <- colnames(meta_data)
  available_screens <- grep(
    "sc[a-zA-Z]+$|DEGAS$",
    colnames(meta_data),
    value = TRUE
  )
  if (!group_by %chin% all_screen_types) {
    Abort(
      "Grouping variable not found in metadata.",
      tips = "Current: {.val {group_by}}\nAvailable grouping variables: {.val {all_screen_types}}"
    )
  }
  # Check available screen types in the Seurat object
  if (is.null(screen_type)) {
    pattern_detect <- paste0(paste(names(ScreenStrategy), collapse = "$|"), "$")
    screen_type <- grep(
      pattern_detect,
      all_screen_types,
      value = TRUE,
      ignore.case = TRUE
    )
  }
  if (!all(screen_type %chin% all_screen_types)) {
    # We don't use `available_screens` for the sake of compatibility to other groups
    Abort(
      "Screen type(s) not found in metadata.",
      tips = "Available screen types: {.val {available_screens}}"
    )
  }

  plot_color <- plot_color %||%
    c(
      "Neutral" = "#CECECE",
      "Other" = "#CECECE",
      "Positive" = "#ff3333",
      "Negative" = "#386c9b"
    )

  # Function to create plot for single screen type
  SinglePlot <- function(single_screen_type, title_suffix = "") {
    stats_df <- meta_data |>
      dplyr::count(
        !!dplyr::sym(group_by),
        !!dplyr::sym(single_screen_type)
      ) |>
      complete_counts(
        !!dplyr::sym(group_by),
        !!dplyr::sym(single_screen_type)
      ) |>
      dplyr::group_by(!!dplyr::sym(group_by)) |>
      dplyr::mutate(Total = sum(`n`)) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        Fraction = ifelse(Total == 0, 0, `n` / Total)
      )

    plot_order <- stats_df |>
      dplyr::filter(!!dplyr::sym(single_screen_type) == "Positive") |>
      dplyr::arrange(dplyr::desc(Fraction)) |>
      dplyr::pull(!!dplyr::sym(group_by))

    # Get label type from misc
    label_type <- screened_seurat@misc[[grep(
      glue::glue("{single_screen_type}_type"),
      names(screened_seurat@misc),
      value = TRUE
    )[1]]]

    if (is.null(label_type)) {
      label_type <- single_screen_type
    }

    # Filter null records
    plot_df <- if (!show_null) {
      dplyr::filter(stats_df, Fraction > 0)
    } else {
      stats_df
    }

    # Create plot title
    current_title <- if (
      length(plot_title) > 1 && length(plot_title) == length(screen_type)
    ) {
      plot_title[which(screen_type == single_screen_type)]
    } else if (title_suffix != "") {
      glue::glue("{plot_title} - {title_suffix}")
    } else {
      plot_title
    }

    plot <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = factor(!!dplyr::sym(group_by), levels = plot_order),
        y = `Fraction`,
        fill = !!dplyr::sym(single_screen_type)
      )
    ) +

      ggplot2::geom_col(position = "stack", width = stack_width) +
      ggplot2::scale_y_continuous(
        labels = function(x) paste0(round(x * 100L, 0L), "%"),
        expand = c(0, 0),
        breaks = seq(0, 1, 0.1)
      ) +
      ggplot2::scale_fill_manual(
        values = plot_color
      ) +
      ggplot2::theme_classic(base_size = 14L) +
      ggplot2::labs(
        x = x_lab,
        y = y_lab,
        fill = label_type,
        title = current_title
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = x_text_angle,
          hjust = 1,
          vjust = 1
        ),
        axis.text = ggplot2::element_text(color = "black"),
        axis.line = ggplot2::element_line(linewidth = axis_linewidth),
        !!!theme_args
      )

    list(stats = stats_df, plot = plot)
  }

  # Process single or multiple screen types
  if (length(screen_type) == 1) {
    # Single screen type - maintain backward compatibility
    result <- SinglePlot(screen_type)

    if (show_plot) {
      methods::show(result$plot)
    }

    return(result)
  }
  # Multiple screen types
  if (verbose) {
    cli::cli_inform(
      "Creating plots for {.val {length(screen_type)}} screen types..."
    )
  }

  # Create plots for all screen types in parallel if possible
  plot_results <- lapply(screen_type, function(st) {
    SinglePlot(st, title_suffix = st)
  })

  # Extract stats and plots
  stats_list <- lapply(plot_results, function(x) x$stats)
  names(stats_list) <- screen_type

  plots_list <- lapply(plot_results, function(x) x$plot)
  names(plots_list) <- screen_type

  # Combine plots using patchwork
  combined_plot <- patchwork::wrap_plots(
    plots_list,
    ncol = ncol,
    nrow = nrow
  )

  if (show_plot) {
    methods::show(combined_plot)
  }

  list(
    stats = stats_list,
    plot = plots_list,
    combined_plot = combined_plot
  )
}

#' @title A Function to Replace `tidyr::complete()`
#' @description
#' Detects if `tidyr::complete()` is available and uses it if so.
#' Otherwise, uses `expand.grid()` and `dplyr::left_join()` to achieve the same result.
#'
#' @keywords internal
complete_counts <- function(data, ..., fill = list(n = 0)) {
  cols <- rlang::enquos(...)
  col_names <- purrr::map_chr(cols, rlang::quo_name)

  if (rlang::is_installed("tidyr")) {
    # usee tidyr::complete
    return(getExportedValue("tidyr", "complete")(
      data,
      !!!cols,
      fill = fill
    ))
  }
  # use expand.grid + left_join
  unique_vals <- purrr::map(col_names, ~ unique(data[[.x]]))

  all_combos <- do.call(
    expand.grid,
    c(
      unique_vals,
      list(stringsAsFactors = FALSE)
    )
  ) |>
    stats::setNames(col_names)

  dplyr::left_join(all_combos, data, by = col_names) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(names(fill)),
        ~ dplyr::coalesce(.x, fill[[dplyr::cur_column()]])
      )
    )
}
