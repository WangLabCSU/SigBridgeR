# IsSkewedDynamic - basic behavior ---------------------------------------

describe("IsSkewedDynamic - basic behavior", {
  it("returns FALSE when observed proportion matches expected_p", {
    # 80% zeros, expected_p = 0.8 -> not skewed
    x <- c(rep(0, 80), rep(1, 20))
    result <- IsSkewedDynamic(x, target = 0, expected_p = 0.8)

    expect_false(result)
  })

  it("returns TRUE when observed proportion deviates significantly", {
    # 50% zeros, expected_p = 0.8 -> difference = 0.3 > n_sd * SD
    x <- c(rep(0, 50), rep(1, 50))
    result <- IsSkewedDynamic(x, target = 0, expected_p = 0.8)

    expect_true(result)
  })

  it("returns NA for empty vector", {
    result <- IsSkewedDynamic(numeric(0))

    expect_true(is.na(result))
  })

  it("returns FALSE for exactly matching proportion (large n)", {
    # 80 zeros out of 100 -> exactly 0.8
    x <- c(rep(0, 80), rep(1, 20))
    result <- IsSkewedDynamic(x, target = 0, expected_p = 0.8, n_sd = 4)

    expect_false(result)
  })

  it("returns TRUE when deviation exceeds n_sd threshold", {
    # n=100, expected_p=0.5, SD = sqrt(0.5*0.5/100) = 0.05
    # With n_sd=2: threshold = 0.1
    # 70% zeros -> diff = 0.2 > 0.1 -> TRUE
    x <- c(rep(0, 70), rep(1, 30))
    result <- IsSkewedDynamic(x, target = 0, expected_p = 0.5, n_sd = 2)

    expect_true(result)
  })

  it("returns FALSE when deviation is within n_sd threshold", {
    # n=100, expected_p=0.5, SD = 0.05
    # With n_sd=3: threshold = 0.15
    # 60% zeros -> diff = 0.1 < 0.15 -> FALSE
    x <- c(rep(0, 60), rep(1, 40))
    result <- IsSkewedDynamic(x, target = 0, expected_p = 0.5, n_sd = 3)

    expect_false(result)
  })

  it("higher n_sd makes it harder to be skewed", {
    # n=100, expected_p=0.5, SD=0.05
    # With n_sd=2: threshold=0.1, 61% zeros -> diff=0.11 > 0.1 -> TRUE
    # With n_sd=5: threshold=0.25, diff=0.11 < 0.25 -> FALSE
    x <- c(rep(0, 61), rep(1, 39))

    result_2sd <- IsSkewedDynamic(x, target = 0, expected_p = 0.5, n_sd = 2)
    result_5sd <- IsSkewedDynamic(x, target = 0, expected_p = 0.5, n_sd = 5)

    expect_true(result_2sd)
    expect_false(result_5sd)
  })
})

# IsSkewedDynamic - input types ------------------------------------------

describe("IsSkewedDynamic - input types", {
  it("works with numeric vector", {
    x <- c(rep(0, 80), rep(1, 20))
    result <- IsSkewedDynamic(x)

    expect_type(result, "logical")
    expect_length(result, 1)
  })

  it("works with integer vector", {
    x <- c(rep(0L, 80), rep(1L, 20))
    result <- IsSkewedDynamic(x)

    expect_type(result, "logical")
    expect_false(result)
  })

  it("works with logical vector", {
    x <- c(rep(FALSE, 80), rep(TRUE, 20))
    result <- IsSkewedDynamic(x, target = 0)

    expect_type(result, "logical")
    expect_false(result)
  })

  it("works with logical vector, target = 1 (TRUE)", {
    x <- c(rep(TRUE, 80), rep(FALSE, 20))
    result <- IsSkewedDynamic(x, target = 1, expected_p = 0.8)

    # 80% TRUE, expected_p=0.8 -> not skewed
    expect_false(result)
  })

  it("treats logical TRUE as 1 and FALSE as 0", {
    x <- c(rep(FALSE, 50), rep(TRUE, 50))
    result <- IsSkewedDynamic(x, target = 0, expected_p = 0.8)

    # 50% FALSE(0), expected_p=0.8 -> skewed
    expect_true(result)
  })
})

# IsSkewedDynamic - target parameter -------------------------------------

describe("IsSkewedDynamic - target parameter", {
  it("uses target = 0 by default", {
    x <- c(rep(0, 50), rep(1, 50))

    result_default <- IsSkewedDynamic(x, expected_p = 0.8)
    result_explicit <- IsSkewedDynamic(x, target = 0, expected_p = 0.8)

    expect_equal(result_default, result_explicit)
  })

  it("counts non-zero target correctly", {
    # 80% ones, target=1, expected_p=0.8 -> not skewed
    x <- c(rep(1, 80), rep(0, 20))
    result <- IsSkewedDynamic(x, target = 1, expected_p = 0.8)

    expect_false(result)
  })

  it("handles target value not present in data", {
    # No 2s in the vector, p_hat=0, expected_p=0.5 -> diff=0.5 > threshold -> TRUE
    x <- c(rep(0, 80), rep(1, 20))
    result <- IsSkewedDynamic(x, target = 2, expected_p = 0.5)

    expect_true(result)
  })

  it("handles fractional target value", {
    x <- c(rep(0.5, 80), rep(1.0, 20))
    result <- IsSkewedDynamic(x, target = 0.5, expected_p = 0.8)

    expect_false(result)
  })
})

# IsSkewedDynamic - expected_p parameter ---------------------------------

describe("IsSkewedDynamic - expected_p parameter", {
  it("uses expected_p = 0.8 by default", {
    x <- c(rep(0, 80), rep(1, 20))

    result_default <- IsSkewedDynamic(x)
    result_explicit <- IsSkewedDynamic(x, expected_p = 0.8)

    expect_equal(result_default, result_explicit)
  })

  it("expected_p = 0.5 with balanced data is not skewed", {
    x <- c(rep(0, 50), rep(1, 50))
    result <- IsSkewedDynamic(x, target = 0, expected_p = 0.5)

    expect_false(result)
  })

  it("expected_p = 0 is valid (all non-target expected)", {
    # expected_p=0 means we expect no target values
    # 20% target vs expected 0% -> skewed
    x <- c(rep(0, 80), rep(1, 20))
    result <- IsSkewedDynamic(x, target = 1, expected_p = 0, n_sd = 2)

    expect_true(result)
  })

  it("expected_p = 1 is valid (all target expected)", {
    # expected_p=1 means we expect all target values
    # 80% target vs expected 100% -> skewed
    x <- c(rep(1, 80), rep(0, 20))
    result <- IsSkewedDynamic(x, target = 1, expected_p = 1, n_sd = 2)

    expect_true(result)
  })
})

# IsSkewedDynamic - n_sd parameter ---------------------------------------

describe("IsSkewedDynamic - n_sd parameter", {
  it("uses n_sd = 4 by default", {
    x <- c(rep(0, 80), rep(1, 20))

    result_default <- IsSkewedDynamic(x)
    result_explicit <- IsSkewedDynamic(x, n_sd = 4)

    expect_equal(result_default, result_explicit)
  })

  it("n_sd = 0 makes any deviation skewed", {
    # n_sd=0 -> threshold=0 -> any diff > 0 -> TRUE (if not exact)
    x <- c(rep(0, 51), rep(1, 49))
    result <- IsSkewedDynamic(x, target = 0, expected_p = 0.5, n_sd = 0)

    # p_hat=0.51, diff=0.01 > 0 -> TRUE
    expect_true(result)
  })

  it("n_sd = 0 with exact match is not skewed", {
    x <- c(rep(0, 50), rep(1, 50))
    result <- IsSkewedDynamic(x, target = 0, expected_p = 0.5, n_sd = 0)

    # p_hat=0.5, diff=0 -> not > 0 -> FALSE
    expect_false(result)
  })
})

# IsSkewedDynamic - edge cases -------------------------------------------

describe("IsSkewedDynamic - edge cases", {
  it("returns NA for empty numeric vector", {
    expect_true(is.na(IsSkewedDynamic(numeric(0))))
  })

  it("returns NA for empty integer vector", {
    expect_true(is.na(IsSkewedDynamic(integer(0))))
  })

  it("returns NA for empty logical vector", {
    expect_true(is.na(IsSkewedDynamic(logical(0))))
  })

  it("handles single-element vector matching target", {
    x <- 0
    result <- IsSkewedDynamic(x, target = 0, expected_p = 0.8)

    # n=1, p_hat=1, SD=sqrt(0.16)=0.4, threshold=4*0.4=1.6
    # diff=0.2 < 1.6 -> FALSE
    expect_false(result)
  })

  it("handles single-element vector not matching target", {
    x <- 1
    result <- IsSkewedDynamic(x, target = 0, expected_p = 0.8)

    # n=1, p_hat=0, diff=0.8, SD=0.4, threshold=1.6 -> FALSE
    expect_false(result)
  })

  it("handles very large vector", {
    n <- 10000
    x <- c(rep(0, as.integer(n * 0.85)), rep(1, as.integer(n * 0.15)))
    result <- IsSkewedDynamic(x, target = 0, expected_p = 0.8)

    # n=10000, p_hat=0.85, SD=sqrt(0.16/10000)=0.004
    # threshold=4*0.004=0.016, diff=0.05 > 0.016 -> TRUE
    expect_true(result)
  })

  it("handles all target values", {
    x <- rep(0, 100)
    result <- IsSkewedDynamic(x, target = 0, expected_p = 0.8)

    # p_hat=1, diff=0.2, SD=0.04, threshold=0.16 -> TRUE
    expect_true(result)
  })

  it("handles no target values", {
    x <- rep(1, 100)
    result <- IsSkewedDynamic(x, target = 0, expected_p = 0.8)

    # p_hat=0, diff=0.8 > 0.16 -> TRUE
    expect_true(result)
  })

  it("larger n makes small deviations detectable", {
    # Small n: diff same but SD larger -> not skewed
    x_small <- c(rep(0, 8), rep(1, 2))  # n=10, 80% zeros
    result_small <- IsSkewedDynamic(x_small, target = 0, expected_p = 0.8)

    # Large n: same proportion but SD smaller -> may be skewed
    x_large <- c(rep(0, 800), rep(1, 200))  # n=1000, 80% zeros
    result_large <- IsSkewedDynamic(x_large, target = 0, expected_p = 0.8)

    # Both are exactly 0.8, neither should be skewed
    expect_false(result_small)
    expect_false(result_large)
  })
})

# IsSkewedDynamic - NA handling ------------------------------------------

describe("IsSkewedDynamic - NA handling", {
  it("aborts when numeric vector contains NA", {
    x <- c(1, 2, NA, 4)

    expect_error(
      IsSkewedDynamic(x),
      "x contains"
    )
  })

  it("aborts when integer vector contains NA", {
    x <- c(1L, 2L, NA_integer_, 4L)

    expect_error(
      IsSkewedDynamic(x),
      "x contains"
    )
  })

  it("aborts when logical vector contains NA", {
    x <- c(TRUE, FALSE, NA)

    expect_error(
      IsSkewedDynamic(x),
      "x contains"
    )
  })
})

# IsSkewedDynamic - symmetry ---------------------------------------------

describe("IsSkewedDynamic - symmetry", {
  it("symmetric: deviation above and below are both detected", {
    # Above expected: 90% zeros vs expected 80%
    x_above <- c(rep(0, 90), rep(1, 10))
    # Below expected: 70% zeros vs expected 80%
    x_below <- c(rep(0, 70), rep(1, 30))

    result_above <- IsSkewedDynamic(x_above, target = 0, expected_p = 0.8)
    result_below <- IsSkewedDynamic(x_below, target = 0, expected_p = 0.8)

    # Both have same absolute deviation (0.1), both should be skewed
    expect_equal(result_above, result_below)
  })

  it("swapping target and non-target inverts the result symmetrically", {
    # 80% zeros, target=0 -> p_hat=0.8 -> matches expected -> FALSE
    x <- c(rep(0, 80), rep(1, 20))
    r1 <- IsSkewedDynamic(x, target = 0, expected_p = 0.8)

    # 20% ones, target=1 -> p_hat=0.2 -> diff=0.6 from 0.8 -> TRUE
    r2 <- IsSkewedDynamic(x, target = 1, expected_p = 0.8)

    expect_false(r1)
    expect_true(r2)
  })
})
