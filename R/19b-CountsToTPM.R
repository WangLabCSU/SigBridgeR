#' Convert bulk RNA-seq raw counts to TPM
#'
#' @param counts A base matrix or Matrix object. Rows are genes, columns are samples.
#' @param gene_length A named numeric vector or named list. Names are gene IDs,
#'   values are gene lengths in bp.
#'
#' @return TPM matrix. Base matrix input returns base matrix;
#'   Matrix input returns Matrix object.
#' @export
CountsToTPM <- function(counts, gene_length) {
  if (is.list(gene_length)) {
    gene_length <- unlist(gene_length, recursive = TRUE, use.names = TRUE)
  }

  gl_names <- names(gene_length)

  if (is.null(gl_names)) {
    Abort("`gene_length` must be named.")
  }

  if (any(nzchar(names(gene_length)))) {
    Abort("`gene_length` contains {.val NA} or empty names.")
  }

  if (any(!is.finite(gene_length) | gene_length <= 0)) {
    Abort("`gene_length` must contain positive finite gene lengths in bp.")
  }

  if (is.data.frame(counts)) {
    counts <- as.matrix(counts)
  }

  if (!inherits(counts, c("matrix", "dgCMatrix", "dgeMatrix"))) {
    cls_counts <- class(counts)
    Abort(
      "Expected `counts` to be a {.cls matrix/dgCMatrix/dgeMatrix}",
      "Current input is {.cls {cls_counts}}"
    )
  }
  CountsToTPM_impl(counts, gene_length)
}
