#' Check whether a matrix looks like bulk RNA-seq raw counts
#'
#' @description
#' Determines whether a numeric matrix is likely to be a bulk RNA-seq raw
#' counts matrix. A matrix is accepted only if it is non-empty, finite,
#' non-negative, integer-like in a sufficiently high proportion, and has a
#' positive library size for every sample. Useful for validating input before
#' running count-based downstream analyses.
#'
#' @param x Numeric matrix. Rows are genes, columns are samples.
#' @param verbose Logical. Whether to print diagnostic statistics (fractions of
#'   non-negative, integer-like, and zero values, library sizes, and the reason
#'   for rejection) to the console. Default: `TRUE`.
#' @param integer_tol Numeric in `[0, 1]`. Absolute tolerance used to decide
#'   whether a value is integer-like (`abs(x - round(x)) <= integer_tol`).
#'   Default: `1e-8`.
#' @param min_integer_fraction Numeric in `[0, 1]`. Minimum fraction of
#'   integer-like values required for the matrix to be considered counts.
#'   Default: `0.95`.
#'
#' @return A single logical value: `TRUE` if `x` is likely a bulk counts
#'   matrix, `FALSE` otherwise.
#'
#' @details
#' `IsCountsMatrix()` rejects a matrix for any of the following reasons:
#'
#' * the matrix is empty;
#' * the matrix contains `NA`, `NaN`, or infinite values;
#' * the matrix contains negative values;
#' * the fraction of integer-like values is below `min_integer_fraction`;
#' * at least one sample (column) has a library size (column sum) `<= 0`.
#'
#' When `verbose = TRUE`, the reason for rejection (or diagnostic statistics on
#' acceptance) is printed to the console.
#'
#' @examples
#' # A typical bulk RNA-seq raw counts matrix (genes x samples)
#' counts <- matrix(
#'   c(100, 200, 300, 150, 250, 350),
#'   nrow = 3, ncol = 2,
#'   dimnames = list(c("GENE1", "GENE2", "GENE3"), c("Sample1", "Sample2"))
#' )
#' IsCountsMatrix(counts)
#'
#' # A log-normalized matrix is not counts
#' IsCountsMatrix(log1p(counts))
#'
#' # Matrices containing negative values are not counts
#' IsCountsMatrix(counts - 1000)
#'
#' # Inspect diagnostic statistics
#' IsCountsMatrix(counts, verbose = TRUE)
#'
#' @family bulk_preprocessing
#'
#' @export
IsCountsMatrix <- function(
  x,
  verbose = TRUE,
  integer_tol = 1e-8,
  min_integer_fraction = 0.95
) {
  chk::chk_matrix(x, x_name = "RNA expression matrix")
  chk::chk_logical(verbose, x_name = "verbose")
  chk::chk_range(integer_tol, x_name = "integer_tol")
  chk::chk_range(min_integer_fraction, x_name = "min_integer_fraction")
  IsCountsMatrixImpl(
    x = x,
    verbose = verbose,
    integer_tol = integer_tol,
    min_integer_fraction = min_integer_fraction
  )
}
