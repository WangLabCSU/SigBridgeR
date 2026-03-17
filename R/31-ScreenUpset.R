# ---- Screened cell intersections ----

#' @title ScreenUpset - Visualize cell type intersections from screened Seurat object
#'
#' @description
#' This function creates an UpSet plot to visualize intersections between different
#' screening methods (e.g., scissor, scAB, scPAS, scPP) in a Seurat object metadata.
#' It calculates all possible combinations of screening types and displays the
#' number of cells positive for each combination.
#'
#' @param screened_seurat A Seurat object containing screening results in metadata.
#'                        Must contain columns with screening types marked as "Positive".
#' @param screen_type Character vector of screening types to analyze.
#'                    Default: NULL, indicating that a self-search pattern will be used.
#' @param show_plot Whether to show the upset plot. Default: FALSE
#' @param n_intersections Number of intersections to display in the plot. Default: 20.
#' @param x_lab Label for the x-axis. Default: "Screen Set Intersections".
#' @param y_lab Label for the y-axis. Default: "Number of Cells".
#' @param title Plot title. Default: "Cell Counts Across Screen Set Intersections".
#' @param bar_color Color for the bars in the plot. Default: "#4E79A7".
#' @param combmatrix_point_color Color for points in the combination matrix. Default: "black".
#' @param verbose Logical, whether to print a message
#' @param ... Additional arguments passed to `ggplot2::theme()` for customizing the plot appearance.
#'
#' @return A list containing two elements:
#'   - `plot`: A ggplot object displaying the UpSet plot
#'   - `stats`: A data frame with intersection statistics including:
#'     - `intersection`: Name of the intersection
#'     - `sets`: List of screening types in the intersection
#'     - `count`: Number of cells positive for all screening types in the intersection
#'
#' @details
#' The function performs the following steps:
#' 1. Validates input parameters and checks if specified screening types exist in metadata
#' 2. Generates all possible combinations of screening types (from 1 to the total number of types)
#' 3. Creates a logical matrix of positive cells for efficient computation
#' 4. Calculates cell counts for each intersection using vectorized operations
#' 5. Creates an UpSet plot visualization with customizable appearance
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' result <- ScreenUpset(screened_seurat = my_seurat_obj)
#'
#' # Customize screening types and appearance
#' result <- ScreenUpset(
#'   screened_seurat = my_seurat_obj,
#'   screen_type = c("scissor", "scAB"),
#'   n_intersections = 15,
#'   title = "Custom Title",
#'   bar_color = "darkred"
#' )
#' }
#'
#'
#' @export
#' @family visualization_function
#'
ScreenUpset <- function(
  screened_seurat,
  screen_type = NULL,
  show_plot = FALSE,
  n_intersections = 20,
  x_lab = "Screen Set Intersections",
  y_lab = "Number of Cells",
  title = "Cell Counts Across Screen Set Intersections",
  bar_color = "#4E79A7",
  combmatrix_point_color = "black",
  verbose = SigBridgeRUtils::getFuncOption("verbose"),
  ...
) {
  rlang::check_installed(c("tibble", "ggplot2", "ggupset"))
  # Robust
  chk::chk_is(screened_seurat, "Seurat")
  if (!is.null(screen_type)) {
    chk::chk_character(screen_type)
  }
  chk::chk_whole_number(n_intersections)
  chk::chk_flag(show_plot)
  purrr::walk(
    list(x_lab, y_lab, title, bar_color, combmatrix_point_color),
    ~ chk::chk_character
  )

  meta_data <- screened_seurat[[]]
  all_screen_types <- colnames(meta_data)

  if (is.null(screen_type)) {
    pattern_detect <- paste0(paste(names(ScreenStrategy), collapse = "$|"), "$")
    screen_type <- grep(
      pattern_detect,
      all_screen_types,
      value = TRUE,
      ignore.case = TRUE
    )
  }
  if (!all(screen_type %in% all_screen_types)) {
    cli::cli_abort(c(
      "x" = "Screen type(s) not found in metadata.",
      ">" = "{.val {screen_type[!screen_type %in% all_screen_types]}}"
    ))
  }

  max_comb <- length(screen_type)

  # Generate all combinations from 1 to max_comb
  all_combinations <- unlist(
    lapply(seq_len(max_comb), function(k) {
      if (length(screen_type) >= k) {
        combs <- utils::combn(screen_type, k, simplify = FALSE)
        stats::setNames(
          combs,
          vapply(
            X = combs,
            FUN = function(comb) {
              if (length(comb) == 1) {
                comb
              } else {
                paste(comb, collapse = " & ")
              }
            },
            FUN.VALUE = character(1)
          )
        )
      }
    }),
    recursive = FALSE
  )

  # Precompute logical matrix for faster filtering
  positive_matrix <- as.matrix(
    meta_data[, screen_type, drop = FALSE] == "Positive"
  )

  # Calculate intersection counts using vectorized operations
  counts <- vapply(
    X = all_combinations,
    FUN = function(sets) {
      # Find rows where all specified sets are positive using vectorized rowSums
      row_matches <- Matrix::rowSums(positive_matrix[,
        sets,
        drop = FALSE
      ]) ==
        length(sets)

      sum(row_matches, na.rm = TRUE)
    },
    FUN.VALUE = numeric(1)
  )

  # Create result data frame
  intersection_data <- tibble::tibble(
    intersection = names(all_combinations),
    sets = all_combinations, # Use I() to preserve list structure
    count = counts
  )

  # Arguments allocated to ggplot2::theme() and ggupset::theme_combmatrix()
  dots <- rlang::list2(...)
  dots$combmatrix.panel.point.color.fill <- dots$combmatrix.panel.point.color.fill %||%
    combmatrix_point_color
  dots$combmatrix.label.make_space <- dots$combmatrix.label.make_space %||%
    FALSE
  verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption("verbose")

  theme_args <- SigBridgeRUtils::FilterArgs4Func(dots, ggplot2::theme)
  combmatrix_args <- SigBridgeRUtils::FilterArgs4Func(
    dots,
    ggupset::theme_combmatrix
  )

  # Create UpSet plot
  upset <- ggplot2::ggplot(
    intersection_data,
    ggplot2::aes(x = `sets`, y = `count`)
  ) +
    ggplot2::geom_col(fill = bar_color, alpha = 0.9, width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = `count`), vjust = -0.5) +
    ggupset::scale_x_upset(
      order_by = "degree",
      sets = screen_type,
      n_intersections = n_intersections
    ) +
    ggplot2::labs(
      x = x_lab,
      y = y_lab,
      title = title
    ) +
    rlang::exec(ggupset::theme_combmatrix, !!!combmatrix_args) +
    ggplot2::theme_minimal() +
    rlang::exec(ggplot2::theme, !!!theme_args)

  if (show_plot) {
    methods::show(upset)
  }

  if (verbose) {
    cli::cli_alert_success("ScreenUpset completed.")
  }

  list(plot = upset, stats = intersection_data)
}
