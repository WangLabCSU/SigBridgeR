# Tests for FPKMToTPM() (R/12d-FPKMToTPM.R).
#
# Public contract:
#   * FPKM values are rescaled per sample so that each column sums to 1e6.
#   * `data` may be a base numeric/integer/logical matrix, a data frame of
#     numeric columns, or an S4 Matrix object; the values are handed to
#     beachmat::initializeCpp() and read column by column in the C++ backend.
#   * The return value is always a dense base matrix with the same dimensions and
#     dimnames as `data`.
#   * Missing values can be treated as zero (na_as_zero = TRUE) or rejected.
#   * Samples whose total FPKM is zero yield an all-zero column and a warning.
#
# Notes on assertions:
#   * cli_abort()/Abort() renders backticks as inline markup, so error patterns
#     must not contain backticks.
#   * Errors raised by the C++ backend come from beachmat/tatami and are not
#     part of this package's API, so those tests only assert that an error is
#     thrown, without pinning the (version-dependent) message.
#   * `...` is forwarded to CheckNA() only when `data` contains missing values
#     AND `verbose = TRUE`.

ref_fpkm_to_tpm <- function(fpkm, na_as_zero = TRUE) {
  fpkm <- as.matrix(fpkm)
  storage.mode(fpkm) <- "double"

  if (na_as_zero) {
    fpkm[is.na(fpkm)] <- 0L
  }

  totals <- colSums(fpkm)
  out <- fpkm
  for (j in seq_len(ncol(fpkm))) {
    out[, j] <- fpkm[, j] * (1e6 / totals[j])
  }
  unname(out)
}

test_that("FPKMToTPM() rescales every sample to one million", {
  fpkm <- matrix(c(1L, 2L, 3L, 4L), nrow = 2L, ncol = 2L)

  tpm <- FPKMToTPM(fpkm, verbose = FALSE)

  expect_equal(colSums(tpm), rep(1e6, 2L))
})

test_that("FPKMToTPM() returns a dense base matrix with the input dimensions", {
  fpkm <- matrix(
    c(1L, 2L, 3L, 4L, 5L, 6L),
    nrow = 3L,
    ncol = 2L,
    dimnames = list(c("g1", "g2", "g3"), c("s1", "s2"))
  )

  tpm <- FPKMToTPM(fpkm, verbose = FALSE)

  expect_true(is.matrix(tpm))
  expect_false(isS4(tpm))
  expect_equal(dim(tpm), dim(fpkm))
  # The C++ layer returns a bare NumericMatrix; dimnames are restored on the R
  # side, matching CountsToTPM().
  expect_equal(dimnames(tpm), dimnames(fpkm))
})

test_that("FPKMToTPM() returns unnamed output when the input is unnamed", {
  fpkm <- matrix(c(1L, 2L, 3L, 4L), nrow = 2L, ncol = 2L)

  expect_null(dimnames(FPKMToTPM(fpkm, verbose = FALSE)))
})

test_that("FPKMToTPM() computes the expected value for each gene", {
  fpkm <- matrix(c(1L, 3L, 6L, 2L), nrow = 2L, ncol = 2L)

  tpm <- FPKMToTPM(fpkm, verbose = FALSE)

  # Column 1 sums to 4, column 2 sums to 8.
  expect_equal(tpm[1L, ], c(1L / 4L * 1e6, 6L / 8L * 1e6))
  expect_equal(tpm[2L, ], c(3L / 4L * 1e6, 2L / 8L * 1e6))
})

test_that("FPKMToTPM() matches a base R reference implementation", {
  set.seed(42L)
  fpkm <- matrix(runif(400L), nrow = 100L, ncol = 4L)

  tpm <- FPKMToTPM(fpkm, verbose = FALSE)

  expect_equal(tpm, ref_fpkm_to_tpm(fpkm), tolerance = 1e-12L)
})

test_that("FPKMToTPM() handles single gene, single sample and 1x1 inputs", {
  one_gene <- matrix(c(2L, 4L, 6L), nrow = 1L, ncol = 3L)
  expect_equal(
    FPKMToTPM(one_gene, verbose = FALSE),
    matrix(1e6, nrow = 1L, ncol = 3L)
  )

  one_sample <- matrix(c(1L, 1L, 2L, 2L), nrow = 4L, ncol = 1L)
  expect_equal(
    FPKMToTPM(one_sample, verbose = FALSE),
    ref_fpkm_to_tpm(one_sample),
    tolerance = 1e-12L
  )

  scalar <- matrix(7L, nrow = 1L, ncol = 1L)
  expect_equal(
    FPKMToTPM(scalar, verbose = FALSE),
    matrix(1e6, nrow = 1L, ncol = 1L)
  )
})

# ---- accepted input types ---------------------------------------------------

test_that("FPKMToTPM() accepts integer and logical matrices", {
  int_fpkm <- matrix(c(1L, 3L, 6L, 2L), nrow = 2L, ncol = 2L)
  expect_equal(
    FPKMToTPM(int_fpkm, verbose = FALSE),
    ref_fpkm_to_tpm(int_fpkm),
    tolerance = 1e-12L
  )

  log_fpkm <- matrix(c(TRUE, FALSE, TRUE, TRUE), nrow = 2L, ncol = 2L)
  expect_equal(
    FPKMToTPM(log_fpkm, verbose = FALSE),
    ref_fpkm_to_tpm(log_fpkm),
    tolerance = 1e-12L
  )
})

test_that("FPKMToTPM() accepts a data frame of numeric columns", {
  df <- data.frame(sample_1 = c(1L, 3L), sample_2 = c(6L, 2L))

  tpm <- FPKMToTPM(df, verbose = FALSE)

  expect_true(is.matrix(tpm))
  expect_equal(dim(tpm), c(2L, 2L))
  expect_equal(dimnames(tpm), dimnames(as.matrix(df)))
  expect_equal(unname(tpm), ref_fpkm_to_tpm(df), tolerance = 1e-12L)
})

test_that("FPKMToTPM() accepts S4 Matrix objects and returns a base matrix", {
  skip_if_not_installed("Matrix")

  values <- matrix(c(1L, 0L, 0L, 2L, 3L, 0L), nrow = 2L, ncol = 3L)
  sparse <- Matrix::Matrix(values, sparse = TRUE)
  dense <- Matrix::Matrix(values, sparse = FALSE)
  triplet <- methods::as(sparse, "TsparseMatrix")

  expect_s4_class(sparse, "dgCMatrix")
  expect_s4_class(dense, "dgeMatrix")

  expected <- ref_fpkm_to_tpm(values)

  # S4 Matrix objects report dimnames as a list of NULLs even when unnamed, so
  # values are compared without dimnames here.
  tpm_sparse <- FPKMToTPM(sparse, verbose = FALSE)
  expect_true(is.matrix(tpm_sparse))
  expect_equal(unname(tpm_sparse), expected, tolerance = 1e-12L)

  expect_equal(
    unname(FPKMToTPM(dense, verbose = FALSE)),
    expected,
    tolerance = 1e-12L
  )

  # Triplet matrices go through beachmat's "unknown matrix" fallback, which
  # emits an informational message before delegating to as.matrix().
  expect_equal(
    unname(suppressMessages(FPKMToTPM(triplet, verbose = FALSE))),
    expected,
    tolerance = 1e-12L
  )
})

test_that("FPKMToTPM() preserves the dimnames of a named S4 Matrix", {
  skip_if_not_installed("Matrix")

  values <- matrix(
    c(1L, 0L, 0L, 2L, 3L, 0L),
    nrow = 2L,
    ncol = 3L,
    dimnames = list(c("g1", "g2"), c("s1", "s2", "s3"))
  )
  sparse <- Matrix::Matrix(values, sparse = TRUE)

  tpm <- FPKMToTPM(sparse, verbose = FALSE)

  expect_equal(dimnames(tpm), dimnames(values))
  expect_equal(unname(tpm), ref_fpkm_to_tpm(values), tolerance = 1e-12L)
})

# ---- missing values ---------------------------------------------------------

test_that("FPKMToTPM() treats NA as zero when na_as_zero = TRUE", {
  fpkm <- matrix(c(1L, NA, 3L, 4L), nrow = 2L, ncol = 2L)

  tpm <- FPKMToTPM(fpkm, na_as_zero = TRUE, verbose = FALSE)

  # Column 1 is (1, NA) so it sums to 1; column 2 is (3, 4) so it sums to 7.
  expect_equal(tpm[1L, ], c(1e6, 3L / 7L * 1e6))
  expect_equal(tpm[2L, ], c(0L, 4L / 7L * 1e6))
  expect_false(anyNA(tpm))
  expect_true(all(is.finite(tpm)))
})

test_that("FPKMToTPM() rejects NA when na_as_zero = FALSE", {
  fpkm <- matrix(c(1L, NA, 3L, 4L), nrow = 2L, ncol = 2L)

  expect_error(
    FPKMToTPM(fpkm, na_as_zero = FALSE, verbose = FALSE),
    "missing value"
  )
  expect_error(
    FPKMToTPM(fpkm, na_as_zero = FALSE, verbose = TRUE),
    "na_as_zero"
  )
})

test_that("FPKMToTPM() treats NaN as missing", {
  fpkm <- matrix(c(1L, NaN, 3L, 4L), nrow = 2L, ncol = 2L)

  tpm <- FPKMToTPM(fpkm, na_as_zero = TRUE, verbose = FALSE)

  expect_equal(tpm, ref_fpkm_to_tpm(fpkm), tolerance = 1e-12L)
})

test_that("FPKMToTPM() sets an all-NA column to zero and warns", {
  fpkm <- matrix(c(NA_real_, NA_real_, 3L, 4L), nrow = 2L, ncol = 2L)

  expect_warning(tpm <- FPKMToTPM(fpkm, verbose = FALSE), "total FPKM of zero")
  expect_equal(tpm[, 1L], c(0L, 0L))
  expect_equal(tpm[, 2L], c(3L / 7L * 1e6, 4L / 7L * 1e6))
})

# ---- zero-sum samples -------------------------------------------------------

test_that("FPKMToTPM() sets a zero-total sample to zero and warns", {
  fpkm <- matrix(c(0L, 0L, 3L, 4L), nrow = 2L, ncol = 2L)

  expect_warning(tpm <- FPKMToTPM(fpkm, verbose = FALSE), "total FPKM of zero")

  expect_equal(tpm[, 1L], c(0L, 0L))
  expect_equal(tpm[, 2L], c(3L / 7L * 1e6, 4L / 7L * 1e6))
})

test_that("FPKMToTPM() warns once for an all-zero matrix", {
  fpkm <- matrix(0L, nrow = 3L, ncol = 2L)

  expect_warning(tpm <- FPKMToTPM(fpkm, verbose = FALSE), "total FPKM of zero")

  expect_equal(tpm, matrix(0L, nrow = 3L, ncol = 2L))
})

# ---- input validation -------------------------------------------------------

test_that("FPKMToTPM() rejects inputs that are not 2-D", {
  # The is_2d() guard rejects anything without a 2-D `dim` attribute before the
  # backend is reached, so these all fail with the package's own error.
  expect_error(
    FPKMToTPM(NULL),
    "2D matrix-like object"
  )
  expect_error(FPKMToTPM(1L:4L), "2D matrix-like object")
  expect_error(FPKMToTPM(list(a = 1L, b = 2L)), "2D matrix-like object")
  expect_error(FPKMToTPM(1L:4L, verbose = FALSE), "2D matrix-like object")
})

test_that("FPKMToTPM() rejects matrix-like input the backend cannot read", {
  # A character matrix is 2-D, so it passes the guard and reaches the backend,
  # which only accepts integer/real storage.
  expect_error(FPKMToTPM(matrix(letters[1L:4L], nrow = 2L), verbose = FALSE))
})

test_that("FPKMToTPM() rejects non-finite values", {
  pos_inf <- matrix(c(1L, Inf, 3L, 4L), nrow = 2L, ncol = 2L)
  expect_error(FPKMToTPM(pos_inf, verbose = FALSE), "finite")

  neg_inf <- matrix(c(1L, -Inf, 3L, 4L), nrow = 2L, ncol = 2L)
  expect_error(FPKMToTPM(neg_inf, verbose = FALSE), "finite")
})

test_that("FPKMToTPM() rejects negative FPKM values", {
  fpkm <- matrix(c(1L, -1L, 3L, 4L), nrow = 2L, ncol = 2L)

  expect_error(FPKMToTPM(fpkm, verbose = FALSE), "negative")
})

test_that("FPKMToTPM() validates na_as_zero and verbose as flags", {
  fpkm <- matrix(c(1L, 2L, 3L, 4L), nrow = 2L, ncol = 2L)

  expect_error(FPKMToTPM(fpkm, na_as_zero = "yes"), class = "chk_error")
  expect_error(FPKMToTPM(fpkm, na_as_zero = 1L), class = "chk_error")
  expect_error(FPKMToTPM(fpkm, na_as_zero = NULL), class = "chk_error")
  expect_error(
    FPKMToTPM(fpkm, na_as_zero = c(TRUE, FALSE)),
    class = "chk_error"
  )

  expect_error(FPKMToTPM(fpkm, verbose = "yes"), class = "chk_error")
  expect_error(FPKMToTPM(fpkm, verbose = 1L), class = "chk_error")
  expect_error(FPKMToTPM(fpkm, verbose = NULL), class = "chk_error")
})

# ---- messaging and `...` forwarding -----------------------------------------

test_that("FPKMToTPM() reports progress when verbose = TRUE", {
  fpkm <- matrix(c(1L, 2L, 3L, 4L), nrow = 2L, ncol = 2L)

  expect_message(FPKMToTPM(fpkm), "TPM conversion completed")
  expect_message(FPKMToTPM(fpkm), "2 genes x 2 samples")
})

test_that("FPKMToTPM() is silent when verbose = FALSE", {
  fpkm <- matrix(c(1L, 2L, 3L, 4L), nrow = 2L, ncol = 2L)

  expect_silent(FPKMToTPM(fpkm, verbose = FALSE))
})

test_that("FPKMToTPM() reports the number of replaced missing values", {
  fpkm <- matrix(c(1L, NA, 3L, 4L), nrow = 2L, ncol = 2L)

  expect_message(
    FPKMToTPM(fpkm, na_as_zero = TRUE, verbose = TRUE),
    "Replacing 1 missing value"
  )
})

test_that("FPKMToTPM() forwards `...` to CheckNA() only for NA input", {
  fpkm <- matrix(c(1L, NA, 3L, 4L), nrow = 2L, ncol = 2L)

  # max_print must be a non-negative integer; CheckNA() aborts on -1, which
  # proves that `...` reached it.
  expect_error(
    FPKMToTPM(fpkm, na_as_zero = TRUE, verbose = TRUE, max_print = -1L),
    class = "chk_error"
  )

  # No missing values -> CheckNA() is never called, so the bad `...` is ignored.
  clean <- matrix(c(1L, 2L, 3L, 4L), nrow = 2L, ncol = 2L)
  expect_silent(FPKMToTPM(clean, verbose = FALSE, max_print = -1L))
  expect_message(FPKMToTPM(clean, verbose = TRUE, max_print = -1L), "completed")
})
