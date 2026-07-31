#' @title Aggregate Rows or Columns with Duplicate Names
#'
#' @description
#' These functions collapse duplicated row names (e.g., gene symbols) or column names (e.g., sample IDs)
#' in matrix-like objects by aggregating values using configurable methods. They support:
#' \describe{
#'   \item{Rows}{\code{\link{AggregateDupRows}}: merges rows sharing the same row name.}
#'   \item{Columns}{\code{\link{AggregateDupCols}}: merges columns sharing the same column name.}
#'   \item{Both}{\code{\link{AggregateDups}}: convenience wrapper applying row-then-column aggregation.}
#' }
#' Designed for expression matrices, count tables, or any numeric data where feature/sample duplication occurs.
#' Handles \code{matrix}, \code{data.frame}, and S4 \code{Matrix} classes (e.g. \code{dgCMatrix}) robustly.
#'
#' @param x A numeric matrix-like object (see Details).
#' @param method Character scalar. Aggregation method (see Methods below).
#' @param row_method Aggregation method for rows (in \code{AggregateDups}). If \code{NULL}, inherits \code{method}.
#' @param col_method Aggregation method for columns (in \code{AggregateDups}). If \code{NULL}, inherits \code{method}.
#' @param verbose Whether to print messages
#' @param ... No usage
#'
#' @section Methods:
#' Supported methods (applied column-wise for rows, row-wise for columns):
#' \describe{
#'   \item{\code{"max"}}{Maximum value per group (default).}
#'   \item{\code{"sum"}}{Sum of values per group.}
#'   \item{\code{"mean"}}{Arithmetic mean (uses \code{na.rm = TRUE}).}
#'   \item{\code{"median"}}{Median value.}
#'   \item{\code{"first"}}{First occurrence in original order.}
#' }
#'
#' @section Input Types and Return Types:
#' \tabular{ll}{
#' Input class          \tab Output class (unless noted) \cr
#' \code{matrix}        \tab \code{matrix} \cr
#' \code{data.frame}    \tab \code{data.frame} \cr
#' S4 \code{Matrix}     \tab \code{matrix} (dense) — S4 attributes dropped for generality \cr
#' }
#'
#' Row/column order in output follows *first occurrence* of each unique name in \code{rownames(x)} / \code{colnames(x)}.
#'
#' @return
#' An aggregated object of the same effective type as \code{x}, with unique row/column names.
#'
#' @family duplicate aggregation
#' @rdname aggregate-dups
#' @name aggregate-dups
NULL


#' @rdname aggregate-dups
#' @export
AggregateDupRows <- function(
  x,
  method = c("max", "sum", "mean", "median", "first"),
  verbose = TRUE,
  ...
) {
  method <- arg_match(method)

  row_names <- rownames(x)
  if (is.null(row_names)) {
    Abort("Input must have row names")
  }

  if (!anyDuplicated(row_names) > 0L) {
    if (verbose) {
      cli::cli_alert_info("No duplicated row names found.")
    }
    return(x)
  }

  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  was_df <- is.data.frame(x)

  mat <- to_numeric_matrix(x)

  row_names <- rownames(mat)
  col_names <- colnames(mat)

  res <- aggregate_dup_rows_cpp(
    x = mat,
    row_names = row_names,
    method = method
  )

  dimnames(res) <- list(unique(row_names), col_names)

  if (was_df) {
    res <- as.data.frame(res, check.names = FALSE)
  }

  res
}


#' @rdname aggregate-dups
#' @export
AggregateDupCols <- function(
  x,
  method = c("max", "sum", "mean", "median", "first"),
  verbose = TRUE,
  ...
) {
  method <- arg_match(method)

  col_names <- colnames(x)

  if (is.null(col_names)) {
    Abort("Input must have column names.")
  }

  if (!anyDuplicated(col_names)) {
    if (isTRUE(verbose)) {
      cli::cli_alert_info("No duplicated column names found.")
    }
    return(x)
  }

  was_df <- is.data.frame(x)

  mat <- to_numeric_matrix(x)

  row_names <- rownames(mat)
  col_names <- colnames(mat)

  res <- aggregate_dup_cols_cpp(
    x = mat,
    col_names = col_names,
    method = method
  )

  dimnames(res) <- list(row_names, unique(col_names))

  if (was_df) {
    res <- as.data.frame(res, check.names = FALSE)
  }

  res
}


#' @rdname aggregate-dups
#' @description
#' Convenience wrapper that first aggregates duplicated rows, then duplicated columns.
#' Useful for cleaning matrices where both feature and sample duplication may occur.
#'
#' @param row_method Aggregation method for rows. Defaults to \code{method}.
#' @param col_method Aggregation method for columns. Defaults to \code{method}.
#'
#' @examples
#' # Full deduplication in one step
#' mat <- matrix(1:16, nrow = 4,
#'               dimnames = list(c("TP53", "TP53", "BRCA1", "ACTB"),
#'                             c("S1", "S1", "S2", "S3")))
#' mat
#'
#' AggregateDups(mat, method = "sum")
#'
#' @export
AggregateDups <- function(
  x,
  method = c("max", "sum", "mean", "median", "first"),
  row_method = method,
  col_method = method,
  verbose = TRUE,
  ...
) {
  row_method <- arg_match(
    row_method,
    c("max", "sum", "mean", "median", "first")
  )
  col_method <- arg_match(
    row_method,
    c("max", "sum", "mean", "median", "first")
  )
  x <- AggregateDupCols(
    x = x,
    method = col_method,
    verbose = verbose,
    ...
  )

  AggregateDupRows(
    x = x,
    method = row_method,
    verbose = verbose,
    ...
  )
}


to_numeric_matrix <- function(x) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  } else if (inherits(x, "Matrix")) {
    x <- as.matrix(x)
  } else if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  if (!is.numeric(x)) {
    stop(
      "Input must be numeric or coercible to a numeric matrix.",
      call. = FALSE
    )
  }

  if (!is.double(x)) {
    storage.mode(x) <- "double"
  }

  x
}
