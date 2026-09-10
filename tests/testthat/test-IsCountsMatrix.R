# Tests for IsCountsMatrix() (R/13d-IsCountsMatrix.R).
#
# Public contract:
#   * `x` must be a 2-D matrix-like object (base matrix, data frame, or S4
#     Matrix); anything else is rejected by the is_2d() guard before any other
#     check runs.
#   * A matrix is accepted only if it is non-empty, finite, non-negative,
#     integer-like in at least `min_integer_fraction` of its values, and every
#     sample (column) has a library size (column sum) > 0.
#   * `verbose`, `integer_tol` and `min_integer_fraction` are validated with
#     chk, which signals a "chk_error" condition.
#   * With `verbose = TRUE` the verdict and the diagnostic fractions are
#     reported through cli_alert_info()/cli_alert_success(), i.e. as message
#     conditions.
#
# Notes on assertions:
#   * Abort() prefixes the message with an ANSI-coloured "[INPUT ERROR]" tag, so
#     only the plain-text part of the message is matched.
#   * The values are read through a beachmat pointer. Base numeric/integer/
#     logical matrices, dgCMatrix and dgeMatrix have dedicated readers, whereas
#     data frames and triplet matrices go through beachmat's "unknown matrix"
#     fallback, which emits a message of its own; those calls are wrapped in
#     suppressMessages().
#   * Each call emits at most one message, so message expectations are not
#     ambiguous.

# Reference implementation mirroring the C++ checks, so expected verdicts come
# from the documented contract instead of hand-computed constants.
ref_is_counts_matrix <- function(
  x,
  integer_tol = 1e-8,
  min_integer_fraction = 0.95
) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"

  if (length(x) == 0L) {
    return(FALSE)
  }
  if (any(!is.finite(x))) {
    return(FALSE)
  }
  if (any(x < 0)) {
    return(FALSE)
  }
  if (mean(abs(x - round(x)) <= integer_tol) < min_integer_fraction) {
    return(FALSE)
  }

  all(colSums(x) > 0)
}

# ---- accepted input ---------------------------------------------------------

describe("IsCountsMatrix - accepted input", {
  it("accepts a typical bulk counts matrix", {
    counts <- matrix(
      c(100, 200, 300, 150, 250, 350),
      nrow = 3,
      ncol = 2,
      dimnames = list(c("G1", "G2", "G3"), c("S1", "S2"))
    )

    expect_true(IsCountsMatrix(counts, verbose = FALSE))
  })

  it("returns a single unnamed logical value", {
    counts <- matrix(1:4, nrow = 2, ncol = 2)

    result <- IsCountsMatrix(counts, verbose = FALSE)

    expect_type(result, "logical")
    expect_length(result, 1L)
    expect_null(names(result))
  })

  it("accepts an integer matrix", {
    counts <- matrix(1:4, nrow = 2, ncol = 2)

    expect_type(counts, "integer")
    expect_true(IsCountsMatrix(counts, verbose = FALSE))
  })

  it("accepts a double matrix holding whole numbers", {
    counts <- matrix(c(10, 20, 30, 40), nrow = 2, ncol = 2)

    expect_type(counts, "double")
    expect_true(IsCountsMatrix(counts, verbose = FALSE))
  })

  it("accepts a logical matrix", {
    # A 0/1 matrix is still a valid count matrix.
    counts <- matrix(c(TRUE, FALSE, TRUE, TRUE), nrow = 2, ncol = 2)

    expect_true(IsCountsMatrix(counts, verbose = FALSE))
  })

  it("accepts a single gene or a single sample", {
    expect_true(IsCountsMatrix(matrix(c(5, 7), nrow = 1), verbose = FALSE))
    expect_true(IsCountsMatrix(matrix(c(5, 7, 11), ncol = 1), verbose = FALSE))
  })

  it("accepts a matrix with a small fraction of non-integer values", {
    # 4 non-integer values out of 100 -> integer fraction 0.96 >= 0.95.
    counts <- matrix(rep(1, 100), nrow = 10, ncol = 10)
    counts[1:4] <- c(0.5, 1.5, 2.5, 3.5)

    expect_true(IsCountsMatrix(counts, verbose = FALSE))
  })

  it("agrees with a base R reference implementation", {
    withr::local_seed(42)
    counts <- matrix(rpois(300, lambda = 50), nrow = 20, ncol = 15)

    expect_true(IsCountsMatrix(counts, verbose = FALSE))
    expect_equal(
      IsCountsMatrix(counts, verbose = FALSE),
      ref_is_counts_matrix(counts)
    )
  })
})

# ---- Matrix and data.frame input --------------------------------------------

describe("IsCountsMatrix - Matrix and data.frame input", {
  it("accepts a dgCMatrix", {
    skip_if_not_installed("Matrix")
    counts <- Matrix::Matrix(
      c(100, 0, 200, 0, 300, 0),
      nrow = 3,
      ncol = 2,
      sparse = TRUE
    )

    expect_s4_class(counts, "dgCMatrix")
    expect_true(IsCountsMatrix(counts, verbose = FALSE))
  })

  it("accepts a dgeMatrix", {
    skip_if_not_installed("Matrix")
    counts <- Matrix::Matrix(
      c(100, 200, 300, 150, 250, 350),
      nrow = 3,
      ncol = 2,
      sparse = FALSE
    )

    expect_s4_class(counts, "dgeMatrix")
    expect_true(IsCountsMatrix(counts, verbose = FALSE))
  })

  it("accepts a triplet Matrix through the beachmat fallback", {
    skip_if_not_installed("Matrix")
    counts <- methods::as(
      Matrix::Matrix(
        c(100, 0, 200, 0, 300, 0),
        nrow = 3,
        ncol = 2,
        sparse = TRUE
      ),
      "TsparseMatrix"
    )

    expect_true(suppressMessages(IsCountsMatrix(counts, verbose = FALSE)))
  })

  it("accepts a data.frame of numeric columns", {
    counts <- data.frame(S1 = c(1, 2), S2 = c(3, 4))

    expect_true(suppressMessages(IsCountsMatrix(counts, verbose = FALSE)))
  })
})

# ---- rejection reasons ------------------------------------------------------

describe("IsCountsMatrix - rejection reasons", {
  it("rejects an empty matrix", {
    expect_false(
      IsCountsMatrix(matrix(numeric(0), nrow = 0, ncol = 0), verbose = FALSE)
    )
    expect_false(
      IsCountsMatrix(matrix(numeric(0), nrow = 3, ncol = 0), verbose = FALSE)
    )
  })

  it("rejects a matrix containing NA", {
    counts <- matrix(c(1, NA, 3, 4), nrow = 2, ncol = 2)

    expect_false(IsCountsMatrix(counts, verbose = FALSE))
  })

  it("rejects a matrix containing NaN", {
    counts <- matrix(c(1, NaN, 3, 4), nrow = 2, ncol = 2)

    expect_false(IsCountsMatrix(counts, verbose = FALSE))
  })

  it("rejects a matrix containing infinite values", {
    expect_false(
      IsCountsMatrix(matrix(c(1, Inf, 3, 4), nrow = 2), verbose = FALSE)
    )
    expect_false(
      IsCountsMatrix(matrix(c(1, -Inf, 3, 4), nrow = 2), verbose = FALSE)
    )
  })

  it("rejects a matrix containing negative values", {
    # Column 1 still sums to 4, so negativity is the only violated rule.
    counts <- matrix(c(5, -1, 3, 4), nrow = 2, ncol = 2)

    expect_false(IsCountsMatrix(counts, verbose = FALSE))
  })

  it("rejects a matrix whose integer fraction is too low", {
    # Every column sums to a positive value, so only the integer-like rule is
    # violated.
    counts <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2, ncol = 2)

    expect_false(IsCountsMatrix(counts, verbose = FALSE))
  })

  it("rejects a matrix with a non-positive library size", {
    # Column 1 is all zeros.
    expect_false(
      IsCountsMatrix(matrix(c(0, 0, 3, 4), nrow = 2, ncol = 2), verbose = FALSE)
    )
    # Every column is all zeros.
    expect_false(IsCountsMatrix(matrix(0, nrow = 2, ncol = 2), verbose = FALSE))
  })
})

# ---- tolerance and threshold arguments --------------------------------------

describe("IsCountsMatrix - integer_tol and min_integer_fraction", {
  it("uses integer_tol to decide which values are integer-like", {
    counts <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2, ncol = 2)

    expect_false(IsCountsMatrix(counts, verbose = FALSE))
    # Every value is within 0.5 of an integer.
    expect_true(IsCountsMatrix(counts, verbose = FALSE, integer_tol = 0.5))
  })

  it("accepts any non-negative matrix when min_integer_fraction is 0", {
    counts <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2, ncol = 2)

    expect_true(
      IsCountsMatrix(counts, verbose = FALSE, min_integer_fraction = 0)
    )
  })

  it("requires every value to be integer-like when min_integer_fraction is 1", {
    counts <- matrix(c(1, 2.5, 3, 4), nrow = 2, ncol = 2)

    expect_false(
      IsCountsMatrix(counts, verbose = FALSE, min_integer_fraction = 1)
    )
  })

  it("accepts a threshold exactly equal to the integer fraction", {
    # 3 of 4 values are integer-like -> fraction 0.75.
    counts <- matrix(c(1, 2, 0.5, 4), nrow = 2, ncol = 2)

    expect_true(
      IsCountsMatrix(counts, verbose = FALSE, min_integer_fraction = 0.75)
    )
  })
})

# ---- argument validation ----------------------------------------------------

describe("IsCountsMatrix - argument validation", {
  it("aborts when x is not a 2-D matrix-like object", {
    expect_error(IsCountsMatrix(NULL), "x must be a 2d matrix")
    expect_error(IsCountsMatrix(1:4), "x must be a 2d matrix")
    expect_error(IsCountsMatrix("not a matrix"), "x must be a 2d matrix")
    expect_error(IsCountsMatrix(list(a = 1)), "x must be a 2d matrix")
    expect_error(
      IsCountsMatrix(array(1:8, dim = c(2, 2, 2))),
      "x must be a 2d matrix"
    )
  })

  it("checks the shape of x before the remaining arguments", {
    expect_error(IsCountsMatrix(1:4, verbose = "yes"), "x must be a 2d matrix")
  })

  it("aborts when verbose is not a flag", {
    counts <- matrix(1:4, nrow = 2, ncol = 2)

    expect_error(IsCountsMatrix(counts, verbose = "yes"), class = "chk_error")
    expect_error(IsCountsMatrix(counts, verbose = 1), class = "chk_error")
    expect_error(IsCountsMatrix(counts, verbose = NULL), class = "chk_error")
    expect_error(
      IsCountsMatrix(counts, verbose = c(TRUE, FALSE)),
      class = "chk_error"
    )
  })

  it("aborts when integer_tol is outside [0, 1]", {
    counts <- matrix(1:4, nrow = 2, ncol = 2)

    expect_error(IsCountsMatrix(counts, integer_tol = -0.1), class = "chk_error")
    expect_error(IsCountsMatrix(counts, integer_tol = 1.1), class = "chk_error")
  })

  it("aborts when min_integer_fraction is outside [0, 1]", {
    counts <- matrix(1:4, nrow = 2, ncol = 2)

    expect_error(
      IsCountsMatrix(counts, min_integer_fraction = -0.5),
      class = "chk_error"
    )
    expect_error(
      IsCountsMatrix(counts, min_integer_fraction = 1.5),
      class = "chk_error"
    )
  })

  it("rejects a non-numeric tolerance", {
    counts <- matrix(1:4, nrow = 2, ncol = 2)

    # chk_range() may reject the type itself or let it fall through to the
    # backend, so the condition class is deliberately not pinned here.
    expect_error(IsCountsMatrix(counts, integer_tol = "a"))
  })

  it("rejects a missing tolerance in the backend", {
    counts <- matrix(1:4, nrow = 2, ncol = 2)

    # NA_real_ slips through chk_range(), which ignores missing values, and is
    # caught by the C++ guards instead.
    expect_error(
      IsCountsMatrix(counts, integer_tol = NA_real_),
      "finite non-negative"
    )
    expect_error(
      IsCountsMatrix(counts, min_integer_fraction = NA_real_),
      "between 0 and 1"
    )
  })
})

# ---- verbose reporting ------------------------------------------------------

describe("IsCountsMatrix - verbose reporting", {
  it("reports the diagnostic fractions when the matrix is accepted", {
    counts <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)

    expect_message(IsCountsMatrix(counts), "likely a counts matrix")
  })

  it("reports the reason for each rejection", {
    expect_message(
      IsCountsMatrix(matrix(numeric(0), nrow = 0, ncol = 0)),
      "matrix is empty"
    )
    expect_message(
      IsCountsMatrix(matrix(c(1, NA, 3, 4), nrow = 2, ncol = 2)),
      "NA, NaN, or Inf"
    )
    expect_message(
      IsCountsMatrix(matrix(c(5, -1, 3, 4), nrow = 2, ncol = 2)),
      "contains negative values"
    )
    expect_message(
      IsCountsMatrix(matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2, ncol = 2)),
      "insufficient integer-like values"
    )
    expect_message(
      IsCountsMatrix(matrix(c(0, 0, 3, 4), nrow = 2, ncol = 2)),
      "library size"
    )
  })

  it("stays silent when verbose = FALSE", {
    counts <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)

    expect_no_message(IsCountsMatrix(counts, verbose = FALSE))
  })
})
