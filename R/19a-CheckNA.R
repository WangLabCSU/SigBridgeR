#' @title Check for Missing Values (NA) in Data
#'
#' @description
#' Identifies and reports the presence of missing values (NA) in vectors or
#' two-dimensional data structures. Provides detailed positional information
#' including row/column indices and names when available. Returns structured
#' information for programmatic use while displaying user-friendly messages.
#'
#' @param data A vector, matrix, data frame, or Matrix object to check for NA values.
#' @param max_print Integer. Maximum number of NA positions to display in console output.
#'   Set to \code{0} to suppress position details. Default: \code{5L}.
#' @param ... Additional arguments (currently unused, reserved for future extensions).
#'
#' @return Invisible list containing NA information:
#'   \describe{
#'     \item{\code{count}}{Total number of NA values found}
#'     \item{\code{positions}}{For vectors: integer vector of positions; For 2D data:
#'           data frame with row/column indices and names (if available)}
#'     \item{\code{names}}{For named vectors: character vector of names corresponding to NA positions}
#'   }
#'   Returns empty list components when no NA values are found.
#'
#' @section Output Behavior:
#' \itemize{
#'   \item \strong{No NAs found}: Prints success message and returns empty list
#'   \item \strong{NAs found}: Prints warning with count and up to \code{max_print} positions
#'   \item \strong{Position details}: Includes row/column names when available for easier debugging
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Vector with NAs
#' vec <- c(10, 20, NA, 40, NA, 60)
#' result <- CheckNA(vec)
#' # Output: Found 2 NA values
#' #         Positions: 3, 5
#'
#' # Example 2: Named vector
#' named_vec <- c(Sample1 = 10, Sample2 = NA, Sample3 = 30, Sample4 = NA)
#' result <- CheckNA(named_vec)
#' # Output: Found 2 NA values
#' #         Positions: 2, 4
#' #         with names: Sample2, Sample4
#'
#' # Example 3: Data frame with row and column names
#' df <- data.frame(
#'   Gene1 = c(1.0, 2.0, NA),
#'   Gene2 = c(NA, 2.5, 3.0),
#'   Gene3 = c(1.5, NA, 3.5),
#'   row.names = c("Cell1", "Cell2", "Cell3")
#' )
#' result <- CheckNA(df)
#' # Output: Found 3 NA values in data
#' #         First 3 positions:
#' #         Row 3 (Cell3), Col 1 (Gene1)
#' #         Row 1 (Cell1), Col 2 (Gene2)
#' #         Row 2 (Cell2), Col 3 (Gene3)
#'
#' # Example 4: Matrix without names
#' mat <- matrix(c(1, NA, 3, 4, 5, NA), nrow = 2, ncol = 3)
#' result <- CheckNA(mat)
#' # Output: Found 2 NA values in data
#' #         First 2 positions:
#' #         Row 2, Col 1
#' #         Row 1, Col 3
#'
#' # Example 5: No NAs found
#' clean_data <- c(1, 2, 3, 4, 5)
#' result <- CheckNA(clean_data)
#' # Output: No NA values found in the vector
#'
#' # Example 6: Programmatic use - extract NA count
#' df <- data.frame(a = c(1, NA, 3), b = c(NA, 2, 3))
#' na_info <- CheckNA(df)
#' na_info$count  # Returns: 2
#' na_info$positions  # Returns: data frame with positions
#' }
#'
#' @seealso \code{\link{is.na}}, \code{\link{complete.cases}}, \code{\link{na.omit}}
CheckNA <- function(data, max_print = 5L, ...) {
  chk::chk_integer(max_print)
  # Convert to data.table if it's a 2D structure but not already a data.table
  is_2d_data <- is_2d(data)

  # Initialize output
  na_info <- list()

  # Handle 1D data (vectors)
  if (!is_2d_data) {
    na_logi <- is.na(data)
    na_count <- sum(na_logi)
    na_positions <- which(na_logi)
    na_info$positions <- na_positions
    na_info$count <- na_count

    if (na_count == 0L) {
      cli::cli_alert_success("No NA values found in the data")
      return(invisible(na_info))
    }

    cli::cli_alert_warning("Found {.pkg {na_count}} NA value{?s}")

    if (max_print > 0L) {
      positions_to_show <- na_positions[
        seq_len(min(max_print, na_count))
      ]
      cli::cli_text("Positions: {.pkg {positions_to_show}}")

      if (na_count > max_print) {
        validate_explain(
          "{.pkg {na_count - max_print}} additional positions not shown",
          .envir = rlang::current_env()
        )
      }
    }

    # For named vectors
    if (!is.null(names(data))) {
      na_names <- names(data)[na_positions]
      na_info$names <- na_names

      if (max_print > 0L) {
        names_to_show <- na_names[seq_len(min(max_print, na_count))]
        validate_explain(
          "with names: {.val {names_to_show}}",
          .envir = rlang::current_env()
        )
      }
    }
    return(invisible(na_info))
  }

  # Handle 2D data (data.frames, matrices, data.tables)

  if (!inherits(data, "Matrix")) {
    data <- Matrix::Matrix(as.matrix(data))
  }

  na_mat <- is.na(data)
  na_count <- sum(na_mat)
  na_info$count <- na_count

  if (na_count == 0L) {
    cli::cli_alert_success("No NA values found in the data")
    return(invisible(na_info))
  }

  cli::cli_alert_warning(
    "Found {.pkg {na_count}} NA value{?s} in data"
  )

  if (max_print > 0L) {
    na_positions <- Matrix::which(na_mat, arr.ind = TRUE) # base matrix
    na_positions <- as.data.frame(na_positions)

    if (!is.null(colnames(data))) {
      na_positions$col_name <- colnames(data)[na_positions$col]
      has_col_names <- TRUE
    } else {
      has_col_names <- FALSE
    }
    if (!is.null(rownames(data))) {
      na_positions$row_name <- rownames(data)[na_positions$row]
      has_row_names <- TRUE
    } else {
      has_row_names <- FALSE
    }
    na_info$positions <- na_positions

    if (max_print > 0L) {
      positions_to_show <- min(max_print, na_count)

      cli::cli_text(
        "First {.pkg {positions_to_show}} position{?s}:"
      )

      if (has_row_names && has_col_names) {
        for (i in seq_len(positions_to_show)) {
          validate_explain(
            "Row {na_positions$row[i]} ({.val {na_positions$row_name[i]}}), \
          Col {na_positions$col[i]} ({.val {na_positions$col_name[i]}})",
            .envir = rlang::current_env()
          )
        }
      } else if (has_row_names) {
        for (i in seq_len(positions_to_show)) {
          validate_explain(
            "Row {na_positions$row[i]} ({.val {na_positions$row_name[i]}}), \
          Col {na_positions$col[i]}",
            .envir = rlang::current_env()
          )
        }
      } else if (has_col_names) {
        for (i in seq_len(positions_to_show)) {
          validate_explain(
            "Row {na_positions$row[i]}, \
            Col {na_positions$col[i]} ({.val {na_positions$col_name[i]}})",
            .envir = rlang::current_env()
          )
        }
      } else {
        for (i in seq_len(positions_to_show)) {
          validate_explain(
            "Row {na_positions$row[i]}, \
            Col {na_positions$col[i]}",
            .envir = rlang::current_env()
          )
        }
      }

      if (na_count > max_print) {
        validate_explain(
          "{.pkg {na_count - max_print}} additional positions not shown",
          .envir = rlang::current_env()
        )
      }
    }
  }

  return(invisible(na_info))
}
