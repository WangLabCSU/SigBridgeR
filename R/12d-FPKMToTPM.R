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
    Abort("`data` must be a 2D matrix-like object", type = "[VALUE ERROR]")
  }

  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }

  # 含非数值列的 data frame 会在这里变成字符矩阵
  if (is.matrix(data) && !is.numeric(data)) {
    Abort(
      "`data` must contain only numeric values",
      type = "[VALUE ERROR]"
    )
  }

  # 逐元素校验全部在 C++ 的单趟扫描中完成
  res <- fpkm_to_tpm_cpp(
    beachmat::initializeCpp(data),
    na_as_zero = na_as_zero
  )

  # NA 报告：保持原行为（verbose 且存在 NA 时展示）
  if (verbose && (res$status == 1L || res$na_count > 0)) {
    CheckNA(data = data, ...)
  }

  if (res$status != 0L) {
    where <- get_cell_position(data, res$row, res$col)
    switch(
      as.character(res$status),
      "1" = Abort(
        "`data` contains missing values",
        "Set `na_as_zero = TRUE` to replace them with zero",
        type = "[VALUE ERROR]"
      ),
      "2" = Abort(
        paste0("`data` contains negative values (first at {where})"),
        type = "[VALUE ERROR]"
      ),
      "3" = Abort(
        paste0("`data` contains infinite values (first at {where})"),
        type = "[VALUE ERROR]"
      ),
      "4" = Abort(
        paste0("The FPKM total overflows ({where})"),
        type = "[VALUE ERROR]"
      ),
      "5" = Abort(
        paste0("The FPKM total is too small to normalize ({where})"),
        type = "[VALUE ERROR]"
      )
    )
  }

  result <- res$tpm
  dimnames(result) <- dimnames(data)

  n_zero <- length(res$zero_cols)
  if (n_zero > 0L) {
    shown <- utils::head(res$zero_cols, 5L)
    cn <- colnames(data)
    lab <- if (is.null(cn)) shown else cn[shown]
    cli::cli_warn(c(
      "{n_zero} sample{?s} ha{?s/ve} a total FPKM of zero.",
      i = "Their TPM values are set to zero (e.g. {.val {lab}})."
    ))
  }

  if (verbose) {
    if (res$na_count > 0) {
      cli::cli_alert_info(
        "Replacing {res$na_count} missing value{?s} with zero"
      )
    }
    cli::cli_alert_success(
      "TPM conversion completed: {nrow(result)} genes x {ncol(result)} samples"
    )
  }

  result
}

get_cell_position <- function(data, i, j) {
  esc <- function(x) gsub("([{}])", "\\1\\1", x)
  rn <- rownames(data)
  cn <- colnames(data)
  paste0(
    "row ",
    if (i > 0L) i else "?",
    if (i > 0L && !is.null(rn)) paste0(" '", esc(rn[i]), "'"),
    ", column ",
    j,
    if (!is.null(cn)) paste0(" '", esc(cn[j]), "'")
  )
}
