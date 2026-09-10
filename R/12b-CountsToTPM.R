#' Convert bulk RNA-seq raw counts to TPM
#'
#' @param counts A base matrix or Matrix object. Rows are genes, columns are samples.
#' @param gene_length A named numeric vector or named list. Names are gene IDs,
#'   values are gene lengths in bp.
#'
#' @return TPM matrix. Base matrix input returns base matrix;
#'   Matrix input returns Matrix object.
#'
#' @examples
#' \dontrun{
#' # Create a counts matrix (3 genes x 2 samples)
#' counts <- matrix(
#'   c(100, 200, 300, 150, 250, 350),
#'   nrow = 3, ncol = 2,
#'   dimnames = list(c("GENE1", "GENE2", "GENE3"), c("Sample1", "Sample2"))
#' )
#'
#' # Gene lengths in bp (must be named and match rownames)
#' gene_length <- c(GENE1 = 1000, GENE2 = 2000, GENE3 = 1500)
#'
#' # Convert counts to TPM
#' tpm <- CountsToTPM(counts, gene_length)
#' colSums(tpm)  # Each column should sum to 1e6
#'
#' # Also works with sparse matrices (dgCMatrix)
#' sparse_counts <- Matrix::Matrix(counts, sparse = TRUE)
#' tpm_sparse <- CountsToTPM(sparse_counts, gene_length)
#' }
#'
#' @export
CountsToTPM <- function(counts, gene_length) {
  row_names <- rownames(counts)
  if (is.null(row_names)) {
    Abort("`counts` must have rownames.", type = "[INPUT ERROR]")
  }

  if (is.list(gene_length)) {
    gene_length <- unlist(gene_length, recursive = TRUE, use.names = TRUE)
  }

  gl_names <- names(gene_length)

  if (is.null(gl_names)) {
    Abort("`gene_length` must be a named vector.", type = "[INPUT ERROR]")
  }

  if (!all(nzchar(names(gene_length)))) {
    Abort(
      "`gene_length` contains {.val NA} or {.val empty names}.",
      type = "[INPUT ERROR]"
    )
  }

  if (any(!is.finite(gene_length) | gene_length <= 0)) {
    Abort(
      "`gene_length` must contain positive finite gene lengths in bp.",
      type = "[INPUT ERROR]"
    )
  }

  if (is.data.frame(counts)) {
    counts <- as.matrix(counts)
  }

  was_s4mat <- if (inherits(counts, "Matrix")) {
    TRUE
  } else {
    FALSE
  }

  if (!is_2d(counts)) {
    cls_counts <- class(counts)
    Abort(
      "Expected `counts` to be a {.cls matrix-like} object",
      "Current input is {.cls {cls_counts}}",
      type = "[TYPE ERROR]"
    )
  }
  result <- counts_to_tpm_cpp(
    beachmat::initializeCpp(counts),
    row_names,
    gene_length
  )
  dimnames(result) <- dimnames(counts)

  if (was_s4mat) {
    result <- Matrix::Matrix(result)
  }

  result
}
