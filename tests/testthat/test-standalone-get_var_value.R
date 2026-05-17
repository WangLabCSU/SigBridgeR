# Source the standalone file so that all .gvv_* helpers and get_var_value
# are available during testing.
source(testthat::test_path("../../R/standalone-get_var_value.R"))

test_that("get_var_value works", {
  f <- function(save_path = "./analysis") {
    save_path_new <- file.path(save_path, "res")
    return(save_path_new)
  }
  f_res <- get_var_value("save_path_new", f) # returns "./analysis/res"

  expect_equal(f_res, "./analysis/res")

  g <- function(a = 1, b = 2) {
    c <- a * 2 + b * 3
    d = c^2
    e <<- d - 1
    e
  }
  g_res <- get_var_value("e", g) # returns 63

  expect_equal(g_res, 63)

  h <- function(a = "A", ...) {
    a <- 1
    return(a)
    a <- 2
    a
  }
  h_res <- get_var_value("a", h) # returns 1

  expect_equal(h_res, 1)
})

test_that("compatible with control flow", {
  f2 <- function(a = 1, b = 2) {
    c <- if (a >= 1) {
      a + b
    } else {
      a * b
    }
  }
  f_res <- get_var_value("c", f2) # returns 3
  expect_equal(f_res, 3)

  g2 <- function(a = 2) {
    for (i in 1:8) {
      a <- a * 2
    }
    a
  }

  g_res <- get_var_value("a", g2) # returns 512
  expect_equal(g_res, 512)

  h2 <- function(a = 2) {
    while (a < 10) {
      a <- a * 2
    }
    a
  }
  h_res <- get_var_value("a", h2) # returns 16
  expect_equal(h_res, 16)

  # repeat with break
  i2 <- function(a = 2) {
    repeat {
      a <- a * 2
      if (a > 20) break
    }
    a
  }
  i_res <- get_var_value("a", i2) # returns 32
  expect_equal(i_res, 32)

  # if with assignments in both branches
  j2 <- function(x = 1) {
    if (x > 0) {
      y <- x * 10
    } else {
      y <- x * (-10)
    }
    y
  }
  expect_equal(get_var_value("y", j2), 10)

  k2 <- function(x = -1) {
    if (x > 0) {
      y <- x * 10
    } else {
      y <- x * (-10)
    }
    y
  }
  expect_equal(get_var_value("y", k2), 10)
})

test_that("compatible with tryCatch", {
  f <- function(x = 1) {
    x <- tryCatch(
      x + 1,
      error = -1
    )
  }

  expect_equal(get_var_value("x", f), 2)
})

test_that("compatible with pkg", {
  # subset-assignment with TRUE condition (branch taken)
  complex_func_true <- function(x = matrix(6:14, 3, byrow = TRUE)) {
    if (x[1, 1] > 5) {
      for (i in 1:3) {
        x[i, i] <- x[i, i] * 2
      }
    }
    stats::dist(x)
  }

  expect_equal(
    get_var_value("x", complex_func_true),
    matrix(c(12, 9, 12, 7, 20, 13, 8, 11, 28), 3, 3)
  )

  # Deterministic function call from stats (without stats:: prefix)
  v <- runif(100L)
  a_func <- function(x = v) {
    y = dplyr::case_when(
      x > 0.5 ~ 1,
      x <= 0.5 ~ 0
    )
    y
  }
  # Reset seed just before get_var_value so the internal eval of runif(100)
  # (triggered during formals resolution) generates the same sequence.

  expect_equal(
    get_var_value("y", a_func),
    dplyr::case_when(
      v > 0.5 ~ 1,
      v <= 0.5 ~ 0
    )
  )
})
