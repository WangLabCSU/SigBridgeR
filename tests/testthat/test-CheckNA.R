# CheckNA -----------------------------------------------------------------

describe("CheckNA - vector input", {
  it("returns empty list for vector with no NAs", {
    clean_vec <- c(1L, 2L, 3L, 4L, 5L)
    result <- CheckNA(clean_vec)

    expect_equal(result$count, 0L)
    expect_equal(result$positions, integer(0L))
    expect_null(result$names)
  })

  it("detects NAs in numeric vector", {
    vec <- c(10L, 20L, NA, 40L, NA, 60L)
    result <- CheckNA(vec)

    expect_equal(result$count, 2L)
    expect_equal(result$positions, c(3L, 5L))
    expect_null(result$names)
  })

  it("detects NAs in integer vector", {
    vec <- c(1L, NA, 3L, NA, 5L)
    result <- CheckNA(vec)

    expect_equal(result$count, 2L)
    expect_equal(result$positions, c(2L, 4L))
  })

  it("detects NAs in character vector", {
    vec <- c("a", NA, "c", NA, "e")
    result <- CheckNA(vec)

    expect_equal(result$count, 2L)
    expect_equal(result$positions, c(2L, 4L))
  })

  it("detects NAs in logical vector", {
    vec <- c(TRUE, NA, FALSE, NA)
    result <- CheckNA(vec)

    expect_equal(result$count, 2L)
    expect_equal(result$positions, c(2L, 4L))
  })

  it("detects NAs in complex vector", {
    vec <- c(1L + 2i, NA_complex_, 3L + 0i)
    result <- CheckNA(vec)

    expect_equal(result$count, 1L)
    expect_equal(result$positions, 2L)
  })

  it("does not treat NULL list elements as NA", {
    # NULL in a list is not the same as NA in R;
    # is_one_na returns FALSE for NULL SEXP elements
    vec <- list(1L, NULL, 3L, NULL)
    result <- CheckNA(vec)

    expect_equal(result$count, 0L)
    expect_equal(result$positions, integer(0L))
  })

  it("includes names for named vector with NAs", {
    named_vec <- c(Sample1 = 10L, Sample2 = NA, Sample3 = 30L, Sample4 = NA)
    result <- CheckNA(named_vec)

    expect_equal(result$count, 2L)
    expect_equal(result$positions, c(2L, 4L))
    expect_equal(result$names, c("Sample2", "Sample4"))
  })

  it("handles all-NA vector", {
    vec <- c(NA_real_, NA_real_, NA_real_)
    result <- CheckNA(vec)

    expect_equal(result$count, 3L)
    expect_equal(result$positions, c(1L, 2L, 3L))
  })

  it("handles empty vector", {
    result <- CheckNA(numeric(0L))

    expect_equal(result$count, 0L)
    expect_equal(result$positions, integer(0L))
  })

  it("returns invisibly", {
    expect_invisible(CheckNA(c(1L, 2L, NA)))
  })
})

describe("CheckNA - vector max_print", {
  it("max_print = 0 suppresses position details", {
    vec <- c(10L, NA, 20L, NA, 30L, NA)
    result <- CheckNA(vec, max_print = 0L)

    expect_equal(result$count, 3L)
    expect_equal(result$positions, c(2L, 4L, 6L))
  })

  it("max_print limits displayed positions but returns all", {
    vec <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
    result <- CheckNA(vec, max_print = 2L)

    expect_equal(result$count, 10L)
    expect_equal(length(result$positions), 10L)
  })

  it("max_print greater than NA count shows all", {
    vec <- c(1L, NA, 3L)
    result <- CheckNA(vec, max_print = 100L)

    expect_equal(result$count, 1L)
    expect_equal(result$positions, 2L)
  })
})

describe("CheckNA - 2D data input", {
  it("detects NAs in data.frame and returns structured positions", {
    df <- data.frame(
      Gene1 = c(1.0, 2.0, NA),
      Gene2 = c(NA, 2.5, 3.0),
      Gene3 = c(1.5, NA, 3.5),
      row.names = c("Cell1", "Cell2", "Cell3")
    )
    result <- CheckNA(df)

    expect_equal(result$count, 3L)
    expect_s3_class(result$positions, "data.frame")
    expect_equal(nrow(result$positions), 3L)
    expect_equal(
      colnames(result$positions),
      c("row", "col", "row_name", "col_name")
    )
  })

  it("detects NAs in data.frame without explicit row/col names", {
    # data.frame always has rownames (auto-generated "1","2",...) and colnames,
    # so row_name and col_name columns are always present for data.frames
    df <- data.frame(a = c(1L, NA, 3L), b = c(NA, 2L, 3L))
    result <- CheckNA(df)

    expect_equal(result$count, 2L)
    expect_equal(nrow(result$positions), 2L)
    expect_equal(
      colnames(result$positions),
      c("row", "col", "row_name", "col_name")
    )
    expect_equal(result$positions$row_name, c("2", "1"))
    expect_equal(result$positions$col_name, c("a", "b"))
  })

  it("detects NAs in matrix", {
    mat <- matrix(c(1L, NA, 3L, 4L, 5L, NA), nrow = 2L, ncol = 3L)
    result <- CheckNA(mat)

    expect_equal(result$count, 2L)
    expect_equal(nrow(result$positions), 2L)
    expect_equal(colnames(result$positions), c("row", "col"))
    expect_false("row_name" %in% colnames(result$positions))
  })

  it("detects NAs in matrix with dimnames", {
    # 2x3 matrix column-major: c(1, NA, 3, 4, 5, NA)
    # NAs at position 1 (0-based) = row 2 col 1, position 5 (0-based) = row 2 col 3
    mat <- matrix(c(1L, NA, 3L, 4L, 5L, NA), nrow = 2L, ncol = 3L)
    rownames(mat) <- c("R1", "R2")
    colnames(mat) <- c("C1", "C2", "C3")
    result <- CheckNA(mat)

    expect_equal(result$count, 2L)
    expect_equal(result$positions$row, c(2L, 2L))
    expect_equal(result$positions$col, c(1L, 3L))
    expect_equal(result$positions$row_name, c("R2", "R2"))
    expect_equal(result$positions$col_name, c("C1", "C3"))
  })

  it("detects NAs in data.frame with only row names", {
    df <- data.frame(a = c(1L, NA), b = c(NA, 2L), row.names = c("R1", "R2"))
    result <- CheckNA(df)

    # data.frame always has colnames, so col_name is present too
    expect_equal(
      colnames(result$positions),
      c("row", "col", "row_name", "col_name")
    )
    expect_equal(result$positions$row_name, c("R2", "R1"))
    expect_equal(result$positions$col_name, c("a", "b"))
  })

  it("detects NAs in data.frame with only column names", {
    df <- data.frame(a = c(1L, NA), b = c(NA, 2L))
    result <- CheckNA(df)

    # data.frame always has auto row names and col names
    expect_equal(
      colnames(result$positions),
      c("row", "col", "row_name", "col_name")
    )
    expect_equal(result$positions$row_name, c("2", "1"))
    expect_equal(result$positions$col_name, c("a", "b"))
  })

  it("returns empty positions for data.frame with no NAs", {
    df <- data.frame(a = c(1L, 2L), b = c(3L, 4L))
    result <- CheckNA(df)

    expect_equal(result$count, 0L)
    expect_equal(nrow(result$positions), 0L)
  })

  it("returns empty positions for matrix with no NAs", {
    mat <- matrix(1L:6L, nrow = 2L, ncol = 3L)
    result <- CheckNA(mat)

    expect_equal(result$count, 0L)
    expect_equal(nrow(result$positions), 0L)
  })

  it("returns invisibly for 2D data", {
    df <- data.frame(a = c(1L, NA), b = c(NA, 2L))
    expect_invisible(CheckNA(df))
  })
})

describe("CheckNA - 2D max_print", {
  it("max_print limits displayed 2D positions but returns all", {
    df <- data.frame(a = c(NA, NA, NA), b = c(NA, NA, NA))
    result <- CheckNA(df, max_print = 2L)

    expect_equal(result$count, 6L)
    expect_equal(nrow(result$positions), 6L)
  })

  it("max_print = 0 suppresses 2D position details", {
    df <- data.frame(a = c(1L, NA), b = c(NA, 2L))
    result <- CheckNA(df, max_print = 0L)

    expect_equal(result$count, 2L)
    expect_equal(nrow(result$positions), 2L)
  })
})

describe("CheckNA - input validation", {
  it("aborts when max_print is not integer", {
    expect_error(CheckNA(c(1L, NA), max_print = "5"), class = "chk_error")
  })

  it("aborts when max_print is negative", {
    expect_error(CheckNA(c(1L, NA), max_print = -1L), class = "chk_error")
  })
})

# scan_na_2d --------------------------------------------------------------

describe("scan_na_2d - data.frame", {
  it("returns count, row, col for data.frame with NAs", {
    df <- data.frame(a = c(1L, NA, 3L), b = c(NA, 2L, 3L))
    res <- scan_na_2d(df)

    expect_equal(res$count, 2L)
    expect_equal(res$row, c(2L, 1L))
    expect_equal(res$col, c(1L, 2L))
  })

  it("returns zero count for data.frame with no NAs", {
    df <- data.frame(a = c(1L, 2L), b = c(3L, 4L))
    res <- scan_na_2d(df)

    expect_equal(res$count, 0L)
    expect_equal(res$row, integer(0L))
    expect_equal(res$col, integer(0L))
  })
})

describe("scan_na_2d - matrix", {
  it("returns count, row, col for matrix with NAs", {
    # 2x3 column-major: c(1, NA, 3, 4, 5, NA)
    # NAs at 0-based index 1 (row=1%2=1+1=2, col=1/2=0+1=1)
    # and index 5 (row=5%2=1+1=2, col=5/2=2+1=3)
    mat <- matrix(c(1L, NA, 3L, 4L, 5L, NA), nrow = 2L, ncol = 3L)
    res <- scan_na_2d(mat)

    expect_equal(res$count, 2L)
    expect_equal(res$row, c(2L, 2L))
    expect_equal(res$col, c(1L, 3L))
  })

  it("returns zero count for matrix with no NAs", {
    mat <- matrix(1L:6L, nrow = 2L)
    res <- scan_na_2d(mat)

    expect_equal(res$count, 0L)
    expect_equal(res$row, integer(0L))
    expect_equal(res$col, integer(0L))
  })
})

describe("scan_na_2d - sparseMatrix", {
  it("handles dgCMatrix with NAs in x slot", {
    skip_if_not_installed("Matrix")
    m <- Matrix::Matrix(c(1L, NA, 0L, 2L, NA, 0L), nrow = 2L, sparse = TRUE)
    res <- scan_na_2d(m)

    expect_equal(res$count, 2L)
    expect_equal(length(res$row), 2L)
    expect_equal(length(res$col), 2L)
  })

  it("handles dgCMatrix with no NAs", {
    skip_if_not_installed("Matrix")
    m <- Matrix::Matrix(c(1L, 2L, 0L, 3L, 4L, 0L), nrow = 2L, sparse = TRUE)
    res <- scan_na_2d(m)

    expect_equal(res$count, 0L)
    expect_equal(res$row, integer(0L))
    expect_equal(res$col, integer(0L))
  })

  it("handles dgTMatrix (triplet) with NAs", {
    skip_if_not_installed("Matrix")
    m <- Matrix::Matrix(c(1L, NA, 0L, 2L, NA, 0L), nrow = 2L, sparse = TRUE)
    # Convert to triplet form
    m_triplet <- as(m, "TsparseMatrix")
    res <- scan_na_2d(m_triplet)

    expect_equal(res$count, 2L)
    expect_equal(length(res$row), 2L)
    expect_equal(length(res$col), 2L)
  })

  it("handles lgCMatrix (pattern-only sparse, no x slot)", {
    skip_if_not_installed("Matrix")
    m <- Matrix::Matrix(c(TRUE, FALSE, TRUE, FALSE), nrow = 2L, sparse = TRUE)
    res <- scan_na_2d(m)

    expect_equal(res$count, 0L)
    expect_equal(res$row, integer(0L))
    expect_equal(res$col, integer(0L))
  })

  it("falls back to as.matrix for dense Matrix with unrecognized subclass", {
    skip_if_not_installed("Matrix")
    # Create a dense Matrix and verify it goes through the fallback path
    m <- Matrix::Matrix(c(1L, NA, 3L, 4L), nrow = 2L, sparse = FALSE)
    # Verify it's not a sparseMatrix
    expect_false(inherits(m, "sparseMatrix"))
    # Verify it's a Matrix
    expect_true(inherits(m, "Matrix"))

    res <- scan_na_2d(m)

    expect_equal(res$count, 1L)
    # 2x2 column-major: c(1, NA, 3, 4), NA at 0-based index 1
    # row = 1 % 2 + 1 = 2, col = 1 / 2 + 1 = 1
    expect_equal(res$row, 2L)
    expect_equal(res$col, 1L)
  })
})

describe("scan_na_2d - fallback", {
  it("handles 2D array via as.matrix fallback", {
    # 2x3 column-major: NAs both in row 2
    arr <- array(c(1L, NA, 3L, 4L, 5L, NA), dim = c(2L, 3L))
    res <- scan_na_2d(arr)

    expect_equal(res$count, 2L)
    expect_equal(res$row, c(2L, 2L))
    expect_equal(res$col, c(1L, 3L))
  })
})
