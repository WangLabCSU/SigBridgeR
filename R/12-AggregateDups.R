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
  method <- SigBridgeRUtils::MatchArg(
    method,
    c("max", "sum", "mean", "median", "first")
  )

  row_names <- rownames(x)
  if (is.null(row_names)) {
    cli::cli_abort(c("x" = "Input must have row names"))
  }

  if (!any(duplicated(row_names))) {
    if (verbose) {
      cli::cli_alert_success("No duplicated row names found.")
    }
    return(x)
  }

  is_matrix <- is.matrix(x) && !is.data.frame(x)
  is_dataframe <- is.data.frame(x)
  is_s4matrix <- inherits(x, "Matrix")

  if (is_s4matrix) {
    x <- as.matrix(x)
    is_matrix <- TRUE
  }

  # Use data.table only if data.frame path needed; for pure matrix, vectorize via tapply/split
  # But your current impl is fine — keep for readability if performance is acceptable

  dt <- data.table::as.data.table(x, keep.rownames = "rname")
  value_cols <- names(dt)[-1]

  res_dt <- switch(
    method,
    max = dt[,
      lapply(.SD, max, na.rm = TRUE),
      by = rname,
      .SDcols = value_cols
    ],
    sum = dt[,
      lapply(.SD, sum, na.rm = TRUE),
      by = rname,
      .SDcols = value_cols
    ],
    mean = dt[,
      lapply(.SD, mean, na.rm = TRUE),
      by = rname,
      .SDcols = value_cols
    ],
    median = dt[,
      lapply(.SD, median, na.rm = TRUE),
      by = rname,
      .SDcols = value_cols
    ],
    first = dt[, .SD[1], by = rname, .SDcols = value_cols]
  )

  uniq_rnames <- unique(row_names)
  data.table::setkey(res_dt, rname)
  res_dt <- res_dt[uniq_rnames]

  if (is_matrix) {
    res <- as.matrix(res_dt[, -1])
    rownames(res) <- res_dt$rname
    res
  } else if (is_dataframe) {
    res <- as.data.frame(res_dt[, -1])
    rownames(res) <- res_dt$rname
    res
  } else {
    res_dt
  }
}


#' @rdname aggregate-dups
#' @export
AggregateDupCols <- function(
  x,
  method = c("max", "sum", "mean", "median", "first"),
  verbose = TRUE,
  ...
) {
  method <- SigBridgeRUtils::MatchArg(
    method,
    c("max", "sum", "mean", "median", "first")
  )

  col_names <- colnames(x)
  if (is.null(col_names)) {
    cli::cli_abort(c("x" = "Input must have column names"))
  }

  if (!any(duplicated(col_names))) {
    if (verbose) {
      cli::cli_alert_success("No duplicated column names found.")
    }
    return(x)
  }

  is_matrix <- is.matrix(x) && !is.data.frame(x)
  is_dataframe <- is.data.frame(x)
  orig_rownames <- rownames(x)

  if (inherits(x, "Matrix")) {
    x <- as.matrix(x)
    is_matrix <- TRUE
  }

  uniq_samples <- unique(col_names)
  col_groups <- split(seq_along(col_names), col_names)[uniq_samples]

  res_list <- lapply(col_groups, function(idx) {
    if (length(idx) == 1L) {
      x[, idx, drop = FALSE]
    } else {
      sub_mat <- x[, idx, drop = FALSE]
      switch(
        method,
        sum = SigBridgeRUtils::rowSums3(sub_mat, na.rm = TRUE),
        mean = SigBridgeRUtils::rowMeans3(sub_mat, na.rm = TRUE),
        max = SigBridgeRUtils::rowMaxs3(sub_mat, na.rm = TRUE),
        median = SigBridgeRUtils::rowMedians3(sub_mat, na.rm = TRUE),
        first = sub_mat[, 1L, drop = FALSE]
      )
    }
  })

  res <- do.call(cbind, res_list)
  dimnames(res) <- list(orig_rownames, uniq_samples)

  if (is_dataframe) as.data.frame(res) else res
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
#' AggregateDups(mat, method = "sum")
#' #>       S1 S2 S3
#' #> TP53   5  7  9
#' #> BRCA1  3  7 11
#' #> ACTB   4  8 12
#'
#' @export
AggregateDups <- function(
  x,
  method = c("max", "sum", "mean", "median", "first"),
  row_method = NULL,
  col_method = NULL,
  verbose = TRUE,
  ...
) {
  method <- SigBridgeRUtils::MatchArg(
    method,
    c("max", "sum", "mean", "median", "first")
  )
  if (is.null(row_method)) {
    row_method <- method
  }
  if (is.null(col_method)) {
    col_method <- method
  }

  x <- AggregateDupRows(x, method = row_method, verbose = verbose, ...)
  AggregateDupCols(x, method = col_method, verbose = verbose, ...)
}
