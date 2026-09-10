# Tests for the duplicate-name aggregation family: AggregateDupRows(),
# AggregateDupCols(), and the AggregateDups() convenience wrapper.
#
# Public contract (see R/13a-AggregateDups.R / man/aggregate-dups.Rd):
#   * Rows / columns are grouped by shared names and collapsed with a method.
#   * Output keeps FIRST occurrence order of each unique name.
#   * matrix  -> matrix, data.frame -> data.frame, S4 Matrix -> dense matrix.
#   * When there is nothing duplicated, x is returned unchanged (a message is
#     emitted when verbose = TRUE).

# ---------------------------------------------------------------------------
# Small reference helpers: reproduce the aggregation in base R so expectations
# mirror the C++ semantics rather than hand-computed constants.
# ---------------------------------------------------------------------------

# Reference for AggregateDupRows(): collapse duplicated rows. Output rows are in
# FIRST-occurrence order of each unique row name (same as the C++ Groups).
ref_agg_rows <- function(mat, fun) {
  rn <- rownames(mat)
  uniq <- rn[!duplicated(rn)]
  res <- do.call(
    rbind,
    lapply(uniq, function(g) {
      rows <- which(rn == g)
      # Aggregate each column independently across the rows of the group.
      # apply() over margin 2 is correct even for a single-row group (it then
      # returns that row's values, i.e. a no-op per column).
      apply(mat[rows, , drop = FALSE], 2, fun)
    })
  )
  matrix(res, ncol = ncol(mat), dimnames = list(uniq, colnames(mat)))
}

# Reference for AggregateDupCols(): collapse duplicated columns. Output columns
# are in FIRST-occurrence order of each unique column name.
ref_agg_cols <- function(mat, fun) {
  cn <- colnames(mat)
  uniq <- cn[!duplicated(cn)]
  cols <- lapply(uniq, function(g) {
    idx <- which(cn == g)
    apply(mat[, idx, drop = FALSE], 1, fun)
  })
  do.call(cbind, cols)
}

# Map a public method string to a base R function matching the C++ aggregator
# (C++ drops NA for max/sum/mean/median, so mirror that here).
method_fun <- function(method) {
  switch(
    method,
    max = function(v) suppressWarnings(max(v, na.rm = TRUE)),
    sum = function(v) sum(v, na.rm = TRUE),
    mean = function(v) mean(v, na.rm = TRUE),
    median = function(v) median(v, na.rm = TRUE),
    first = function(v) v[1],
    stop("Unknown method")
  )
}

# Example matrix with duplicated rows AND duplicated columns.
dup_matrix <- function() {
  m <- matrix(
    c(
      1,
      3,
      5,
      7, # row A1 (col1..col4)
      3,
      9,
      5,
      21, # row A2
      5,
      15,
      9,
      35, # row B1
      7,
      21,
      10,
      49 # row C1
    ),
    nrow = 4,
    byrow = TRUE
  )
  dimnames(m) <- list(
    c("G1", "G1", "G2", "G3"),
    c("S1", "S2", "S1", "S3")
  )
  m
}

# ---------------------------------------------------------------------------
# AggregateDupRows()
# ---------------------------------------------------------------------------

test_that("AggregateDupRows collapses duplicated rows with method='sum'", {
  m <- dup_matrix()
  res <- AggregateDupRows(m, method = "sum")

  # Only unique names remain, in first-occurrence order.
  expect_identical(rownames(res), c("G1", "G2", "G3"))

  # Numerical equality against a base-R reconstruction.
  expect_equal(unname(res), unname(ref_agg_rows(m, method_fun("sum"))))
})

test_that("AggregateDupRows honours every aggregation method", {
  m <- dup_matrix()
  for (method in c("max", "sum", "mean", "median", "first")) {
    res <- AggregateDupRows(m, method = method)
    expect_equal(
      unname(res),
      unname(ref_agg_rows(m, method_fun(method))),
      info = paste("method =", method)
    )
    expect_identical(rownames(res), c("G1", "G2", "G3"))
  }
})

test_that("AggregateDupRows keeps a single-column / single-row edge case", {
  m <- matrix(c(1, 5, 2), nrow = 1, dimnames = list("row", c("S1", "S1", "S2")))
  # No duplicated ROWS -> unchanged input.
  expect_identical(AggregateDupRows(m, verbose = FALSE), m)

  # Rows duplicated but single column.
  m2 <- matrix(c(2, 4, 8), nrow = 3, dimnames = list(c("A", "A", "B"), "col"))
  res <- AggregateDupRows(m2, method = "sum", verbose = FALSE)
  expect_equal(unname(res), matrix(c(6, 8), ncol = 1))
  expect_identical(rownames(res), c("A", "B"))
})

test_that("AggregateDupRows returns a data.frame for a data.frame input", {
  df <- data.frame(v1 = c(1, 3, 5), v2 = c(2, 4, 6))
  # data.frame() refuses duplicate row names, so set them directly.
  attr(df, "row.names") <- c("A", "A", "B")
  res <- AggregateDupRows(df, method = "sum", verbose = FALSE)

  expect_s3_class(res, "data.frame")
  expect_identical(rownames(res), c("A", "B"))
  expect_equal(res$v1, c(4, 5))
  expect_equal(res$v2, c(6, 6))
})

test_that("AggregateDupRows returns input unchanged when no rows are duplicated", {
  m <- matrix(1:6, nrow = 3, dimnames = list(c("A", "B", "C"), c("S1", "S2")))

  expect_message(AggregateDupRows(m), "No duplicated row names")
  # quiet path returns identical object
  expect_identical(AggregateDupRows(m, verbose = FALSE), m)
})

test_that("AggregateDupRows errors on missing row names", {
  m <- matrix(1:4, nrow = 2)
  expect_error(AggregateDupRows(m), "row names")
})

test_that("AggregateDupRows rejects an invalid method", {
  m <- dup_matrix()
  expect_error(AggregateDupRows(m, method = "bogus"), "must be one of")
})

# ---------------------------------------------------------------------------
# AggregateDupCols()
# ---------------------------------------------------------------------------

test_that("AggregateDupCols collapses duplicated columns with method='max'", {
  m <- dup_matrix()
  res <- AggregateDupCols(m, method = "max")

  expect_identical(colnames(res), c("S1", "S2", "S3"))
  expect_equal(unname(res), unname(ref_agg_cols(m, method_fun("max"))))
})

test_that("AggregateDupCols honours every aggregation method", {
  m <- dup_matrix()
  for (method in c("max", "sum", "mean", "median", "first")) {
    res <- AggregateDupCols(m, method = method)
    expect_equal(
      unname(res),
      unname(ref_agg_cols(m, method_fun(method))),
      info = paste("method =", method)
    )
    expect_identical(colnames(res), c("S1", "S2", "S3"))
  }
})

test_that("AggregateDupCols returns a data.frame for a data.frame input", {
  df <- data.frame(
    S1a = c(1, 10),
    S1b = c(3, 30),
    S2 = c(5, 50)
  )
  colnames(df) <- c("S1", "S1", "S2")
  rownames(df) <- c("G1", "G2")

  res <- AggregateDupCols(df, method = "sum", verbose = FALSE)
  expect_s3_class(res, "data.frame")
  expect_identical(colnames(res), c("S1", "S2"))
  expect_equal(res[["S1"]], c(4, 40))
  expect_equal(res[["S2"]], c(5, 50))
})

test_that("AggregateDupCols returns input unchanged when no columns duplicate", {
  m <- matrix(1:4, nrow = 2, dimnames = list(c("A", "B"), c("S1", "S2")))

  expect_message(AggregateDupCols(m), "No duplicated column names")
  expect_identical(AggregateDupCols(m, verbose = FALSE), m)
})

test_that("AggregateDupCols errors on missing column names", {
  m <- matrix(1:4, nrow = 2)
  expect_error(AggregateDupCols(m), "column names")
})

test_that("AggregateDupCols rejects an invalid method", {
  m <- dup_matrix()
  expect_error(AggregateDupCols(m, method = "bogus"), "must be one of")
})

# ---------------------------------------------------------------------------
# AggregateDups(): column aggregation then row aggregation
# ---------------------------------------------------------------------------

test_that("AggregateDups equals col-then-row aggregation", {
  m <- dup_matrix()
  res <- AggregateDups(m, method = "sum")

  expect_equal(
    res,
    AggregateDupRows(
      AggregateDupCols(m, method = "sum"),
      method = "sum",
      verbose = FALSE
    )
  )
  expect_identical(rownames(res), c("G1", "G2", "G3"))
  expect_identical(colnames(res), c("S1", "S2", "S3"))
})

test_that("AggregateDups accepts independent row/col methods", {
  m <- dup_matrix()
  res <- AggregateDups(
    m,
    row_method = "max",
    col_method = "sum",
    verbose = FALSE
  )

  expect_equal(
    res,
    AggregateDupRows(
      AggregateDupCols(m, method = "sum", verbose = FALSE),
      method = "max",
      verbose = FALSE
    )
  )
})

test_that("AggregateDups keeps data.frame type", {
  df <- as.data.frame(dup_matrix())
  # as.data.frame() uniquifies duplicated row names ("G1", "G1.1", ...);
  # restore the original duplicates so the aggregation is actually exercised.
  attr(df, "row.names") <- rownames(dup_matrix())
  res <- AggregateDups(df, method = "sum", verbose = FALSE)
  expect_s3_class(res, "data.frame")
  expect_identical(rownames(res), c("G1", "G2", "G3"))
  expect_identical(colnames(res), c("S1", "S2", "S3"))
})

test_that("AggregateDups handles a fully unique matrix (no-op)", {
  m <- matrix(1:4, nrow = 2, dimnames = list(c("A", "B"), c("S1", "S2")))
  expect_identical(AggregateDups(m, verbose = FALSE), m)
})
