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
#'
#'   * **`count`**: Total number of NA values found
#'   * **`positions`**: For vectors: integer vector of positions; For 2D data:
#'     data frame with row/column indices and names (if available)
#'   * **`names`**: For named vectors: character vector of names corresponding to NA positions
#'
#'   Returns empty list components when no NA values are found.
#'
#' @section Output Behavior:
#' * **No NAs found**: Prints success message and returns empty list
#' * **NAs found**: Prints warning with count and up to `max_print` positions
#' * **Position details**: Includes row/column names when available for easier debugging
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
  chk::chk_gte(max_print, 0L)
  check_installed("methods")

  if (!is_2d(data)) {
    res <- check_na_vector_cpp(data)

    na_count <- res$count
    na_positions <- res$positions

    na_info <- list(
      count = na_count,
      positions = na_positions
    )

    if (na_count == 0) {
      cli::cli_alert_success("No NA values found in the data")
      return(invisible(na_info))
    }

    cli::cli_alert_warning(
      "Found {na_count} NA value{if (na_count == 1) '' else 's'}"
    )

    if (max_print > 0L) {
      n_show <- min(max_print, length(na_positions))
      cli::cli_text(
        "Positions: {paste(na_positions[seq_len(n_show)], collapse = ', ')}"
      )

      if (na_count > n_show) {
        cli::cli_text("{na_count - n_show} additional positions not shown")
      }
    }

    if (!is.null(names(data))) {
      na_names <- names(data)[na_positions]
      na_info$names <- na_names

      if (max_print > 0L) {
        n_show <- min(max_print, length(na_names))
        cli::cli_text(
          "with names: {paste(na_names[seq_len(n_show)], collapse = ', ')}"
        )
      }
    }

    return(invisible(na_info))
  }

  res <- scan_na_2d(data)

  na_count <- res$count

  if (na_count == 0) {
    na_info <- list(
      count = na_count,
      positions = data.frame(row = integer(), col = integer())
    )
    cli::cli_alert_success("No NA values found in the data")
    return(invisible(na_info))
  }

  pos <- data.frame(
    row = res$row,
    col = res$col
  )

  rn <- tryCatch(rownames(data), error = function(e) NULL)
  cn <- tryCatch(colnames(data), error = function(e) NULL)

  has_row_names <- !is.null(rn)
  has_col_names <- !is.null(cn)

  if (has_row_names) {
    pos$row_name <- rn[pos$row]
  }

  if (has_col_names) {
    pos$col_name <- cn[pos$col]
  }

  na_info <- list(
    count = na_count,
    positions = pos
  )

  cli::cli_alert_warning(
    "Found {na_count} NA value{if (na_count == 1) '' else 's'} in data"
  )

  if (max_print > 0L) {
    n_show <- min(max_print, nrow(pos))

    cli::cli_text("First {n_show} position{if (n_show == 1) '' else 's'}:")

    for (i in seq_len(n_show)) {
      r <- pos$row[i]
      c <- pos$col[i]

      if (has_row_names && has_col_names) {
        cli::cli_text(
          "Row {r} ({.val {pos$row_name[i]}}), Col {c} ({.val {pos$col_name[i]}})"
        )
      } else if (has_row_names) {
        cli::cli_text(
          "Row {r} ({.val {pos$row_name[i]}}), Col {c}"
        )
      } else if (has_col_names) {
        cli::cli_text(
          "Row {r}, Col {c} ({.val {pos$col_name[i]}})"
        )
      } else {
        cli::cli_text("Row {r}, Col {c}")
      }
    }

    if (na_count > n_show) {
      cli::cli_text("{na_count - n_show} additional positions not shown")
    }
  }

  invisible(na_info)
}

scan_na_2d <- function(data) {
  if (is.data.frame(data)) {
    return(check_na_dataframe_cpp(data, nrow(data)))
  }

  if (is.matrix(data)) {
    return(check_na_dense2d_cpp(data, dim(data)))
  }

  if (inherits(data, "Matrix")) {
    slots <- methods::slotNames(data)
    d <- methods::slot(data, "Dim")

    if (inherits(data, "sparseMatrix")) {
      if (all(c("p", "i") %in% slots)) {
        if (!("x" %in% slots)) {
          return(list(count = 0, row = integer(), col = integer()))
        }

        return(check_na_sparse_csc_cpp(
          methods::slot(data, "x"),
          methods::slot(data, "i"),
          methods::slot(data, "p"),
          d
        ))
      }

      if (all(c("i", "j", "x") %in% slots)) {
        return(check_na_sparse_triplet_cpp(
          methods::slot(data, "x"),
          methods::slot(data, "i"),
          methods::slot(data, "j"),
          d
        ))
      }
    }

    if (all(c("x", "Dim") %in% slots)) {
      xslot <- methods::slot(data, "x")
      if (length(xslot) == prod(d)) {
        return(check_na_dense2d_cpp(xslot, d))
      }
    }

    # fallback for special Matrix classes
    m <- as.matrix(data)
    return(check_na_dense2d_cpp(m, dim(m)))
  }

  m <- as.matrix(data)
  check_na_dense2d_cpp(m, dim(m))
}
