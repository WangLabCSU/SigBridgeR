# compile_lhs -------------------------------------------------------------

describe("compile_lhs - symbols", {
  it("compiles a simple symbol to assign plan", {
    plan <- compile_lhs(quote(x))

    expect_equal(length(plan), 1)
    expect_equal(plan[[1]]$type, "assign")
    expect_equal(plan[[1]]$name, "x")
    expect_equal(plan[[1]]$path, integer(0))
  })

  it("compiles underscore as ignore plan", {
    plan <- compile_lhs(quote(`_`))

    expect_equal(length(plan), 1)
    expect_equal(plan[[1]]$type, "ignore")
    expect_equal(plan[[1]]$path, integer(0))
  })

  it("compiles dot as ignore plan", {
    plan <- compile_lhs(quote(.))

    expect_equal(length(plan), 1)
    expect_equal(plan[[1]]$type, "ignore")
    expect_equal(plan[[1]]$path, integer(0))
  })
})

describe("compile_lhs - c() calls", {
  it("compiles flat c() to multiple assign plans with paths", {
    plan <- compile_lhs(quote(c(x, y, z)))

    expect_equal(length(plan), 3)
    expect_equal(plan[[1]]$type, "assign")
    expect_equal(plan[[1]]$name, "x")
    expect_equal(plan[[1]]$path, 1L)

    expect_equal(plan[[2]]$type, "assign")
    expect_equal(plan[[2]]$name, "y")
    expect_equal(plan[[2]]$path, 2L)

    expect_equal(plan[[3]]$type, "assign")
    expect_equal(plan[[3]]$name, "z")
    expect_equal(plan[[3]]$path, 3L)
  })

  it("compiles nested c() with hierarchical paths", {
    plan <- compile_lhs(quote(c(x, c(y, z))))

    expect_equal(length(plan), 3)
    expect_equal(plan[[1]]$name, "x")
    expect_equal(plan[[1]]$path, 1L)

    expect_equal(plan[[2]]$name, "y")
    expect_equal(plan[[2]]$path, c(2L, 1L))

    expect_equal(plan[[3]]$name, "z")
    expect_equal(plan[[3]]$path, c(2L, 2L))
  })

  it("compiles deeply nested c() calls", {
    plan <- compile_lhs(quote(c(c(a, b), c(c, d))))

    expect_equal(length(plan), 4)
    expect_equal(plan[[1]]$path, c(1L, 1L))
    expect_equal(plan[[2]]$path, c(1L, 2L))
    expect_equal(plan[[3]]$path, c(2L, 1L))
    expect_equal(plan[[4]]$path, c(2L, 2L))
  })

  it("handles mix of ignore and assign in c()", {
    plan <- compile_lhs(quote(c(x, `_`, y)))

    expect_equal(length(plan), 3)
    expect_equal(plan[[1]]$type, "assign")
    expect_equal(plan[[1]]$name, "x")

    expect_equal(plan[[2]]$type, "ignore")

    expect_equal(plan[[3]]$type, "assign")
    expect_equal(plan[[3]]$name, "y")
  })

  it("handles empty c() call", {
    plan <- compile_lhs(quote(c()))

    expect_equal(length(plan), 0)
  })

  it("handles single element c()", {
    plan <- compile_lhs(quote(c(x)))

    expect_equal(length(plan), 1)
    expect_equal(plan[[1]]$name, "x")
    expect_equal(plan[[1]]$path, 1L)
  })
})

describe("compile_lhs - errors", {
  it("throws error for invalid destructuring target (non-symbol, non-c call)", {
    expect_error(
      compile_lhs(quote(1 + 2)),
      "Invalid destructuring target"
    )
  })

  it("throws error for non-c function call", {
    expect_error(
      compile_lhs(quote(list(x, y))),
      "Invalid destructuring target"
    )
  })
})

# %<-% --------------------------------------------------------------------

describe("%<-% - basic assignment", {
  it("assigns single value", {
    c(x) %<-% list(42)

    expect_equal(x, 42)
  })

  it("assigns multiple values from list", {
    c(x, y, z) %<-% list(10, 20, 30)

    expect_equal(x, 10)
    expect_equal(y, 20)
    expect_equal(z, 30)
  })

  it("assigns from atomic vector", {
    c(a, b, c_val) %<-% c(1, 2, 3)

    expect_equal(a, 1)
    expect_equal(b, 2)
    expect_equal(c_val, 3)
  })

  it("assigns from integer vector", {
    c(i, j) %<-% c(10L, 20L)

    expect_equal(i, 10L)
    expect_equal(j, 20L)
  })

  it("assigns from character vector", {
    c(first, second) %<-% c("hello", "world")

    expect_equal(first, "hello")
    expect_equal(second, "world")
  })

  it("assigns from logical vector", {
    c(t_val, f_val) %<-% c(TRUE, FALSE)

    expect_equal(t_val, TRUE)
    expect_equal(f_val, FALSE)
  })

  it("returns rhs invisibly", {
    rhs <- list(1, 2)
    expect_invisible(c(a, b) %<-% rhs)
  })

  it("preserves rhs value in return", {
    rhs <- list(100, 200)
    result <- c(x_ret, y_ret) %<-% rhs
    expect_equal(result, rhs)
  })
})

describe("%<-% - ignore patterns", {
  it("ignores elements with underscore", {
    c(x, `_`, z) %<-% list(10, 20, 30)

    expect_equal(x, 10)
    expect_equal(z, 30)
    expect_false(exists("_", inherits = FALSE))
  })

  it("ignores elements with dot", {
    c(x, ., z) %<-% list(10, 20, 30)

    expect_equal(x, 10)
    expect_equal(z, 30)
  })

  it("ignores first element", {
    c(`_`, y) %<-% list(10, 20)

    expect_equal(y, 20)
  })

  it("ignores last element", {
    c(x, `_`) %<-% list(10, 20)

    expect_equal(x, 10)
  })

  it("ignores all elements", {
    c(`_`, `_`, `_`) %<-% list(1, 2, 3)

    # should not error; nothing assigned
    expect_true(TRUE)
  })
})

describe("%<-% - nested destructuring", {
  it("assigns nested list elements", {
    c(x, c(y, z)) %<-% list(1, list(2, 3))

    expect_equal(x, 1)
    expect_equal(y, 2)
    expect_equal(z, 3)
  })

  it("assigns deeply nested structures", {
    c(c(a, b), c(c_val, d)) %<-% list(list(1, 2), list(3, 4))

    expect_equal(a, 1)
    expect_equal(b, 2)
    expect_equal(c_val, 3)
    expect_equal(d, 4)
  })

  it("handles ignore in nested destructuring", {
    c(x, c(`_`, z)) %<-% list(10, list(20, 30))

    expect_equal(x, 10)
    expect_equal(z, 30)
  })

  it("handles three-level nesting", {
    c(a, c(b, c(c_val, d))) %<-% list(1, list(2, list(3, 4)))

    expect_equal(a, 1)
    expect_equal(b, 2)
    expect_equal(c_val, 3)
    expect_equal(d, 4)
  })
})

describe("%<-% - edge cases", {
  it("assigns single-element destructure", {
    c(x) %<-% list(99)

    expect_equal(x, 99)
  })

  it("assigns to existing variables (overwrites)", {
    x <- 5
    c(x) %<-% list(10)

    expect_equal(x, 10)
  })

  it("works with named lists", {
    c(x, y) %<-% list(a = 1, b = 2)

    expect_equal(x, 1)
    expect_equal(y, 2)
  })

  it("handles empty destructure gracefully", {
    # c() with no args should be a no-op
    c() %<-% list(1, 2, 3)

    expect_true(TRUE)
  })

  it("handles real-world pattern: c(mat, bulk, pheno) %<-% LoadRefData", {
    mock_data <- list(
      matrix(1:6, nrow = 2),
      c(1.0, 2.0, 3.0),
      data.frame(x = 1:3)
    )
    c(mat_exam, bulk, pheno) %<-% mock_data

    expect_equal(dim(mat_exam), c(2, 3))
    expect_equal(bulk, c(1.0, 2.0, 3.0))
    expect_s3_class(pheno, "data.frame")
  })
})
