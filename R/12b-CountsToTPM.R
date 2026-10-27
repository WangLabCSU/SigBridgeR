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
  if (is.null(row_names) || anyNA(row_names)) {
    Abort("`counts` must have non-NA rownames.", type = "[INPUT ERROR]")
  }
  if (!is_2d(counts)) {
    cls_counts <- class(counts)
    Abort(
      "Expected `counts` to be a {.cls matrix-like} object",
      "Current input is {.cls {cls_counts}}",
      type = "[TYPE ERROR]"
    )
  }
  if (is.list(gene_length)) {
    if (any(lengths(gene_length) != 1L)) {
      Abort(
        "Each element of `gene_length` must be a scalar.",
        type = "[INPUT ERROR]"
      )
    }
    gene_length <- unlist(gene_length, recursive = TRUE, use.names = TRUE)
  }

  gl_names <- names(gene_length)
  if (is.null(gl_names) || anyNA(gl_names) || !all(nzchar(gl_names))) {
    Abort(
      "`gene_length` must be a named vector, without {.val NA} or empty names.",
      type = "[INPUT ERROR]"
    )
  }
  if (any(!is.finite(gene_length) | gene_length <= 0L)) {
    Abort(
      "`gene_length` must contain positive finite gene lengths in bp.",
      type = "[INPUT ERROR]"
    )
  }
  if (anyDuplicated(gl_names)) {
    dup <- unique(gl_names[duplicated(gl_names)])
    Abort(
      "`gene_length` has duplicated names, e.g. {.val {utils::head(dup, 5)}}",
      type = "[INPUT ERROR]"
    )
  }

  idx <- match(row_names, gl_names)
  if (anyNA(idx)) {
    missing_genes <- utils::head(row_names[is.na(idx)], 5)
    n_missing <- sum(is.na(idx))
    Abort(
      "`gene_length` misses {n_missing} gene(s), e.g. {.val {missing_genes}}",
      type = "[INPUT ERROR]"
    )
  }
  inv_length <- 1 / as.double(gene_length)[idx]

  if (is.data.frame(counts)) {
    counts <- as.matrix(counts)
  }

  was_s4mat <- if (inherits(counts, "Matrix")) {
    TRUE
  } else {
    FALSE
  }

  result <- counts_to_tpm_cpp(
    beachmat::initializeCpp(counts),
    inv_length
  )
  dimnames(result) <- dimnames(counts)

  if (was_s4mat) {
    result <- Matrix::Matrix(result)
  }

  result
}
