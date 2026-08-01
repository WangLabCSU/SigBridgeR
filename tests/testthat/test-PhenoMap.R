# PhenoMap - vector input -------------------------------------------------

describe("PhenoMap - vector input", {
  it("transforms numeric vector with two conditions", {
    v <- c(5, 15, 25, 35)
    result <- PhenoMap(v, v > 20 ~ "High", v <= 20 ~ "Low")

    expect_equal(result, c("Low", "Low", "High", "High"))
  })

  it("transforms with three conditions and .default", {
    v <- c(5, 15, 25)
    result <- PhenoMap(v, v > 20 ~ "High", v > 10 ~ "Medium", .default = "Low")

    expect_equal(result, c("Low", "Medium", "High"))
  })

  it("returns numeric output when all values are numeric", {
    v <- c(5, 15, 25)
    result <- PhenoMap(v, v > 10 ~ 1, v <= 10 ~ 0)

    expect_equal(result, c(0, 1, 1))
  })

  it("uses .default when no condition matches", {
    v <- c(5, 15, 25)
    result <- PhenoMap(v, v > 30 ~ "Very High", .default = "Other")

    expect_equal(result, c("Other", "Other", "Other"))
  })

  it("returns NA for unmatched when .default is NA", {
    v <- c(5, 15, 25)
    result <- PhenoMap(v, v > 20 ~ "High")

    expect_equal(result, c(NA, NA, "High"))
  })

  it("handles single-element vector", {
    v <- 42
    result <- PhenoMap(v, v > 10 ~ "Big", .default = "Small")

    expect_equal(result, "Big")
  })

  it("preserves names of named vector", {
    v <- c(a = 5, b = 15, c = 25)
    result <- PhenoMap(v, v > 10 ~ "High", .default = "Low")

    expect_equal(names(result), c("a", "b", "c"))
    expect_equal(result, c(a = "Low", b = "High", c = "High"))
  })

  it("handles character vector input", {
    v <- c("apple", "banana", "cherry")
    result <- PhenoMap(
      v,
      nchar(v) > 5 ~ "long",
      .default = "short"
    )

    expect_equal(result, c("short", "long", "long"))
  })

  it("handles logical vector input", {
    v <- c(TRUE, FALSE, TRUE)
    result <- PhenoMap(v, v == TRUE ~ "yes", .default = "no")

    expect_equal(result, c("yes", "no", "yes"))
  })

  it("handles empty vector", {
    v <- numeric(0)
    result <- PhenoMap(v, v > 0 ~ "positive", .default = "other")

    expect_equal(result, character(0))
  })
})

describe("PhenoMap - vector conditions evaluated sequentially", {
  it("first matching condition wins", {
    v <- c(5, 15, 25, 35)
    # v > 10 matches 15, 25, 35; but v > 20 is checked first
    result <- PhenoMap(v, v > 20 ~ "High", v > 10 ~ "Medium", .default = "Low")

    expect_equal(result, c("Low", "Medium", "High", "High"))
  })

  it("overlapping conditions use first match", {
    v <- c(5, 15, 25)
    result <- PhenoMap(
      v,
      v > 10 ~ ">10",
      v > 20 ~ ">20",
      .default = "<=10"
    )

    expect_equal(result, c("<=10", ">10", ">10"))
  })

  it("all conditions matching yields first-condition result", {
    v <- c(100)
    result <- PhenoMap(
      v,
      v > 0 ~ "first",
      v > 50 ~ "second",
      v > 90 ~ "third"
    )

    expect_equal(result, "first")
  })
})

# PhenoMap - data frame input ---------------------------------------------

describe("PhenoMap - data frame input", {
  it("transforms a column in a data.frame in place", {
    df <- data.frame(age = c(25, 35, 45, 55, 65), name = letters[1:5])
    result <- PhenoMap(df, age < 30 ~ "Young", age < 50 ~ "Middle", .default = "Senior")

    expect_s3_class(result, "data.table")
    expect_equal(result$age, c("Young", "Middle", "Middle", "Senior", "Senior"))
    expect_equal(result$name, letters[1:5])
  })

  it("transforms a column with numeric output", {
    df <- data.frame(value = c(5, 15, 25, 35))
    result <- PhenoMap(df, value > 20 ~ 1, value <= 20 ~ 0)

    expect_equal(result$value, c(0, 0, 1, 1))
  })

  it("uses .default for unmatched rows in data.frame", {
    df <- data.frame(score = c(10, 30, 50))
    result <- PhenoMap(df, score > 40 ~ "High", .default = "Not High")

    expect_equal(result$score, c("Not High", "Not High", "High"))
  })

  it("does not modify other columns", {
    df <- data.frame(x = c(1, 2, 3), y = c("a", "b", "c"), z = c(10, 20, 30))
    result <- PhenoMap(df, x > 1 ~ "Big", .default = "Small")

    expect_equal(result$y, c("a", "b", "c"))
    expect_equal(result$z, c(10, 20, 30))
    expect_equal(result$x, c("Small", "Big", "Big"))
  })

  it("handles data.frame with a single row", {
    df <- data.frame(val = 42)
    result <- PhenoMap(df, val > 10 ~ "High", .default = "Low")

    expect_equal(result$val, "High")
  })

  it("handles data.frame with all rows matching .default", {
    df <- data.frame(val = c(1, 2, 3))
    result <- PhenoMap(df, val > 100 ~ "Huge", .default = "Normal")

    expect_equal(result$val, c("Normal", "Normal", "Normal"))
  })

  it("handles data.frame with all rows matching a condition", {
    df <- data.frame(val = c(10, 20, 30))
    result <- PhenoMap(df, val > 5 ~ "Big")

    expect_equal(result$val, c("Big", "Big", "Big"))
  })

  it("works with .default = NA in data.frame (fcase default)", {
    df <- data.frame(val = c(1, 2, 3))
    result <- PhenoMap(df, val > 2 ~ "Big")

    expect_equal(result$val, c(NA, NA, "Big"))
  })

  it("first condition evaluated uses the correct column", {
    df <- data.frame(a = c(1, 2, 3), b = c(10, 20, 30))
    # conditions reference column "a", so only "a" is transformed
    result <- PhenoMap(df, a > 1 ~ "Big", .default = "Small")

    expect_equal(result$a, c("Small", "Big", "Big"))
    expect_equal(result$b, c(10, 20, 30))
  })
})

# PhenoMap - matrix input -------------------------------------------------

describe("PhenoMap - matrix input", {
  it("transforms a matrix column via data.table conversion", {
    mat <- matrix(c(5, 15, 25, 35), ncol = 1, dimnames = list(NULL, "val"))
    result <- PhenoMap(mat, val > 20 ~ "High", .default = "Low")

    expect_equal(result$val, c("Low", "Low", "High", "High"))
  })

  it("handles multi-column matrix (uses first column from conditions)", {
    mat <- matrix(c(5, 15, 10, 20), ncol = 2, dimnames = list(NULL, c("x", "y")))
    result <- PhenoMap(mat, x > 10 ~ "Big", .default = "Small")

    expect_equal(result$x, c("Small", "Big"))
    expect_equal(result$y, c(10, 20))
  })
})

# PhenoMap - data.table input ---------------------------------------------

describe("PhenoMap - data.table input", {
  it("transforms a column in a data.table", {
    dt <- data.table::data.table(val = c(5, 15, 25))
    result <- PhenoMap(dt, val > 10 ~ "High", .default = "Low")

    expect_s3_class(result, "data.table")
    expect_equal(result$val, c("Low", "High", "High"))
  })

  it("preserves other data.table columns", {
    dt <- data.table::data.table(a = c(1, 2, 3), b = c("x", "y", "z"))
    result <- PhenoMap(dt, a > 1 ~ "Big", .default = "Small")

    expect_equal(result$b, c("x", "y", "z"))
    expect_equal(result$a, c("Small", "Big", "Big"))
  })
})

# PhenoMap - input validation ---------------------------------------------

describe("PhenoMap - input validation", {
  it("aborts when no conditions provided", {
    expect_error(
      PhenoMap(c(1, 2, 3)),
      "Condition is empty"
    )
  })

  it("aborts when conditions are not formulas", {
    expect_error(
      PhenoMap(c(1, 2, 3), "not a formula"),
      "Not all conditions are formula"
    )
  })

  it("aborts with mixed formula and non-formula conditions", {
    expect_error(
      PhenoMap(c(1, 2, 3), x > 1 ~ "yes", "not a formula"),
      "Not all conditions are formula"
    )
  })
})

# is_2d -------------------------------------------------------------------

describe("is_2d", {
  it("returns TRUE for data.frame", {
    expect_true(is_2d(data.frame(x = 1:3)))
  })

  it("returns TRUE for matrix", {
    expect_true(is_2d(matrix(1:6, nrow = 2)))
  })

  it("returns TRUE for data.table", {
    expect_true(is_2d(data.table::data.table(x = 1:3)))
  })

  it("returns FALSE for vector", {
    expect_false(is_2d(c(1, 2, 3)))
  })

  it("returns FALSE for list", {
    expect_false(is_2d(list(1, 2, 3)))
  })

  it("returns FALSE for NULL", {
    expect_false(is_2d(NULL))
  })

  it("returns TRUE for 2D array", {
    expect_true(is_2d(array(1:6, dim = c(2, 3))))
  })

  it("returns FALSE for 1D array", {
    expect_false(is_2d(array(1:6, dim = 6)))
  })

  it("returns FALSE for 3D array", {
    expect_false(is_2d(array(1:8, dim = c(2, 2, 2))))
  })
})
