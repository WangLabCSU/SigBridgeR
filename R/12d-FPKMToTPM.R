#' Convert bulk RNA-seq FPKM to TPM
#'
#' @description
#' Converts a gene-level FPKM matrix to TPM (Transcripts Per Million). Because
#' the per-gene length factor is constant within a sample, the FPKM-to-TPM
#' transformation requires no gene lengths: each column is scaled so that its
#' values sum to \eqn{10^6}.
#'
#' For a given sample, every FPKM value is multiplied by
#' \eqn{10^6 / \sum(\text{FPKM})}. Samples whose FPKM column sums to zero
#' produce a zero column in the result and trigger a warning.
#'
#' Missing values (\code{NA}) are either replaced by zero
#' (\code{na_as_zero = TRUE}, the default) or trigger an error
#' (\code{na_as_zero = FALSE}).
#'
#' @param data A two-dimensional matrix-like object of FPKM values (a base
#'   matrix, a data frame, or an S4 \code{Matrix}). Rows are genes, columns are
#'   samples. Data frames are coerced and must contain only numeric columns.
#'   Values must be non-negative and finite.
#' @param na_as_zero A logical flag. If \code{TRUE} (default), missing values
#'   are replaced by zero before normalization. If \code{FALSE}, an error is
#'   raised when \code{data} contains missing values.
#' @param verbose A logical flag. If \code{TRUE} (default), diagnostic
#'   messages (missing-value report and completion summary) are shown.
#' @param ... Additional arguments passed to [CheckNA()] when \code{data}
#'   contains missing values and \code{verbose} is \code{TRUE} (for example
#'   \code{max_print}).
#'
#' @return A numeric TPM matrix with the same dimensions and \code{dimnames} as
#'   \code{data}. Each column sums to \eqn{10^6} unless the corresponding FPKM
#'   column sums to zero.
#'
#' @examples
#' \dontrun{
#' # Create an FPKM matrix (3 genes x 2 samples)
#' fpkm <- matrix(
#'   c(10, 20, 30, 15, 25, 35),
#'   nrow = 3, ncol = 2,
#'   dimnames = list(c("GENE1", "GENE2", "GENE3"), c("Sample1", "Sample2"))
#' )
#'
#' # Convert FPKM to TPM
#' tpm <- FPKMToTPM(fpkm)
#' colSums(tpm)  # Each column should sum to 1e6
#'
#' # Error if missing values are present and na_as_zero = FALSE
#' fpkm_na <- fpkm
#' fpkm_na[1, 1] <- NA
#' tpm2 <- FPKMToTPM(fpkm_na, na_as_zero = FALSE)
#' }
#'
#' @export
FPKMToTPM <- function(data, na_as_zero = TRUE, verbose = TRUE, ...) {
  chk::chk_flag(na_as_zero)
  chk::chk_flag(verbose)

  if (!is_2d(data)) {
    Abort(
      "`data` must be a 2D matrix-like object",
      type = "[VALUE ERROR]"
    )
  }

  # is.infinite()/is.finite() have no method for lists, so a data frame has to
  # be coerced before the finite-value guard below.
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }

  if (any(is.infinite(data), na.rm = TRUE)) {
    Abort(
      "`data` contains infinite values",
      type = "[VALUE ERROR]"
    )
  }

  if (any(data < 0, na.rm = TRUE)) {
    Abort(
      "`data` contains negative values",
      type = "[VALUE ERROR]"
    )
  }

  if (anyNA(data)) {
    if (verbose) {
      CheckNA(data = data, ...)
    }
    if (!na_as_zero) {
      Abort(
        "`data` contains missing values",
        "Set `na_as_zero = TRUE` to replace them with zero",
        type = "[VALUE ERROR]"
      )
    }
  }

  result <- fpkm_to_tpm_cpp(
    beachmat::initializeCpp(data),
    na_as_zero = na_as_zero,
    verbose = verbose
  )

  dimnames(result) <- dimnames(data)
  return(result)
}
