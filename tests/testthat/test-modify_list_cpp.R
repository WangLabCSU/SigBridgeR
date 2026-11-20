# Tests for src/modify_list.cpp — modify_list() and modify_list2(), exercised
# through their R (Rcpp-exported) interfaces.
#
# Run independently with:
#   devtools::load_all()
#   testthat::test_file("tests/testthat/test-modify_list_cpp.R")

# ---- modify_list(): recursive merge, keep_null = TRUE -----------------------

describe("modify_list - recursive merge", {
  it("replaces a named element", {
    expect_equal(
      modify_list(list(a = 1, b = 2), list(a = 10)),
      list(a = 10, b = 2)
    )
  })

  it("appends new named elements", {
    expect_equal(
      modify_list(list(a = 1), list(c = 3)),
      list(a = 1, c = 3)
    )
  })

  it("recursively merges nested lists", {
    expect_equal(
      modify_list(
        list(a = list(p = 1, q = 2)),
        list(a = list(q = 20, r = 30))
      ),
      list(a = list(p = 1, q = 20, r = 30))
    )
  })

  it("replaces a non-list element with a list", {
    expect_equal(
      modify_list(list(a = 1), list(a = list(new = TRUE))),
      list(a = list(new = TRUE))
    )
  })

  it("replaces a list element with a non-list", {
    expect_equal(
      modify_list(list(a = list(p = 1)), list(a = 5)),
      list(a = 5)
    )
  })

  it("appends unnamed elements with an empty name", {
    result <- modify_list(list(a = 1), list(99))
    expect_identical(unname(result), unname(list(1, 99)))
    expect_identical(names(result), c("a", ""))
  })

  it("keeps NULL as a list element with keep_null = TRUE (default)", {
    result <- modify_list(list(a = 1, b = 2), list(a = NULL))
    expect_length(result, 2)
    expect_null(result$a)
    expect_identical(result$b, 2)
    expect_identical(names(result), c("a", "b"))
  })

  it("stores unnamed NULL with keep_null = TRUE", {
    result <- modify_list(list(a = 1), list(NULL), keep_null = TRUE)
    expect_length(result, 2)
    expect_null(result[[2]])
  })

  it("handles empty inputs", {
    expect_equal(modify_list(list(), list()), list())
    expect_equal(modify_list(list(a = 1), list()), list(a = 1))
  })
})

# ---- modify_list(): keep_null = FALSE ---------------------------------------

describe("modify_list - NULL deletion (keep_null = FALSE)", {
  it("deletes a named element whose value is NULL", {
    expect_equal(
      modify_list(list(a = 1, b = 2), list(a = NULL), keep_null = FALSE),
      list(b = 2)
    )
  })

  it("deletes several named elements", {
    expect_equal(
      modify_list(
        list(a = 1, b = 2, c = 3),
        list(a = NULL, c = NULL),
        keep_null = FALSE
      ),
      list(b = 2)
    )
  })

  it("does nothing when the name is absent", {
    expect_equal(
      modify_list(list(a = 1), list(z = NULL), keep_null = FALSE),
      list(a = 1)
    )
  })

  it("ignores unnamed NULL elements", {
    result <- modify_list(list(a = 1), list(NULL), keep_null = FALSE)
    expect_equal(result, list(a = 1))
  })

  it("recursively deletes NULL inside a nested list", {
    expect_equal(
      modify_list(
        list(cfg = list(x = 1, y = 2)),
        list(cfg = list(x = NULL)),
        keep_null = FALSE
      ),
      list(cfg = list(y = 2))
    )
  })
})

# ---- modify_list(): input validation ----------------------------------------

describe("modify_list - input validation", {
  it("errors on a data.frame x (is_plain_list guard)", {
    expect_error(modify_list(data.frame(a = 1), list()), "`x` must be a list")
  })

  it("errors on a data.frame val (is_plain_list guard)", {
    expect_error(modify_list(list(), data.frame(a = 1)), "`val` must be a list")
  })

  it("Rcpp coerces atomic vectors to lists before the guard", {
    result <- modify_list(1:3, list())
    expect_type(result, "list")
    expect_identical(unname(result), list(1L, 2L, 3L))
  })
})

# ---- modify_list2(): val-based merge ----------------------------------------

describe("modify_list2 - val-based merge", {
  it("recursively merges matching named sub-lists", {
    x <- list(a = list(p = 1, q = 2), b = 1, c = list(z = 0))
    val <- list(a = list(q = 20, r = 30), b = list(new = TRUE), d = 4)
    result <- modify_list2(x, val)
    expect_equal(result$a, list(p = 1, q = 20, r = 30))
    expect_equal(result$b, list(new = TRUE))
    expect_equal(result$d, 4)
    expect_identical(names(result), c("a", "b", "d"))
  })

  it("does not carry over x elements absent from val", {
    x <- list(a = 1, b = 2)
    val <- list(a = 10)
    result <- modify_list2(x, val)
    expect_equal(result, list(a = 10))
    expect_length(result, 1)
  })

  it("passes through unnamed val elements", {
    result <- modify_list2(list(1), list(list(a = 1), 2))
    expect_length(result, 2)
    expect_equal(result[[1]], list(a = 1))
    expect_identical(result[[2]], 2)
    expect_null(names(result))
  })

  it("passes through NA-named and empty-named val elements", {
    x <- list(a = list(p = 1), b = list(p = 2))
    val <- setNames(list(list(q = 1), list(q = 2)), c(NA_character_, ""))
    result <- modify_list2(x, val)
    expect_equal(result[[1]], list(q = 1))
    expect_equal(result[[2]], list(q = 2))
  })

  it("uses the first match when x has duplicated names", {
    x <- setNames(list(list(p = 1), list(p = 100)), c("a", "a"))
    val <- list(a = list(q = 2))
    result <- modify_list2(x, val)
    expect_equal(result$a, list(p = 1, q = 2))
  })

  it("replaces when the matching x element is not a list", {
    x <- list(a = 1)
    val <- list(a = list(new = TRUE))
    expect_equal(modify_list2(x, val), list(a = list(new = TRUE)))
  })

  it("keeps only val names on the result", {
    x <- setNames(list(1, 2), c("a", "b"))
    val <- setNames(list(10, 20), c("x", "y"))
    result <- modify_list2(x, val)
    expect_identical(names(result), c("x", "y"))
    expect_identical(unname(result), list(10, 20))
  })

  it("handles empty val", {
    expect_equal(modify_list2(list(a = 1), list()), list())
  })

  it("errors on non-list x", {
    expect_error(modify_list2(1:3, list()), "`x` must be a list")
  })

  it("errors on non-list val", {
    expect_error(modify_list2(list(), 1:3), "`val` must be a list")
  })
})
