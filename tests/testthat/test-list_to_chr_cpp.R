# Tests for src/list_to_chr.cpp — list_to_character_recursive(), exercised
# through its R (Rcpp-exported) interface and the public list_to_chr() wrapper.
#
# Run independently with:
#   devtools::load_all()
#   testthat::test_file("tests/testthat/test-list_to_chr_cpp.R")

describe("list_to_character_recursive - list rendering", {
  it("renders a flat named list", {
    expect_equal(
      list_to_character_recursive(list(a = 1, b = 2)),
      "a: 1, b: 2"
    )
  })

  it("uses 1-based indices for unnamed elements", {
    expect_equal(
      list_to_character_recursive(list(10, 20)),
      "1: 10, 2: 20"
    )
  })

  it("renders nested lists with parentheses", {
    expect_equal(
      list_to_character_recursive(list(a = list(b = 1, c = 2))),
      "a: (b: 1, c: 2)"
    )
  })

  it("handles deeply nested lists", {
    expect_equal(
      list_to_character_recursive(list(a = list(b = list(c = 3)))),
      "a: (b: (c: 3))"
    )
  })

  it("mixes named and unnamed elements", {
    expect_equal(
      list_to_character_recursive(list(a = 1, 2, c = 3)),
      "a: 1, 2: 2, c: 3"
    )
  })

  it("returns an empty string for an empty list", {
    expect_equal(list_to_character_recursive(list()), "")
  })

  it("uses index for empty-string names", {
    x <- setNames(list(1, 2), c("", "b"))
    expect_equal(list_to_character_recursive(x), "1: 1, b: 2")
  })

  it("uses index for NA names", {
    x <- setNames(list(1, 2), c(NA_character_, "b"))
    expect_equal(list_to_character_recursive(x), "1: 1, b: 2")
  })
})

describe("list_to_character_recursive - value types", {
  it("renders NA as NA", {
    expect_equal(list_to_character_recursive(list(x = NA)), "x: NA")
  })

  it("renders NULL as () (Rf_isNewList(NULL) is TRUE on modern R)", {
    expect_equal(list_to_character_recursive(list(x = NULL)), "x: ()")
  })

  it("renders atomic vectors as comma-separated values", {
    expect_equal(
      list_to_character_recursive(list(x = c(1, 2, 3))),
      "x: 1, 2, 3"
    )
  })

  it("renders character vectors containing NA", {
    expect_equal(
      list_to_character_recursive(list(x = c("a", NA, "c"))),
      "x: a, NA, c"
    )
  })

  it("renders logical vectors", {
    expect_equal(
      list_to_character_recursive(list(x = c(TRUE, FALSE, NA))),
      "x: TRUE, FALSE, NA"
    )
  })

  it("renders functions as <function>", {
    expect_equal(
      list_to_character_recursive(list(f = function(x) x)),
      "f: <function>"
    )
  })

  it("renders environments as <environment>", {
    expect_equal(
      list_to_character_recursive(list(e = new.env())),
      "e: <environment>"
    )
  })
})

describe("list_to_chr - public interface", {
  it("returns a single character string matching the C++ wrapper", {
    x <- list(a = 1, b = list(c = 2))
    expect_identical(list_to_chr(x), list_to_character_recursive(x))
    expect_type(list_to_chr(x), "character")
    expect_length(list_to_chr(x), 1)
  })

  it("matches the documented as.character() example", {
    x <- list(
      a = 1,
      b = list(c = 2, d = list(e = "hello")),
      f = c(1, 2, 3)
    )
    expect_equal(
      list_to_chr(x),
      "a: 1, b: (c: 2, d: (e: hello)), f: 1, 2, 3"
    )
  })

  it("renders a list with NULL and NA entries", {
    expect_equal(
      list_to_chr(list("foo", bar = NA, baz = NULL)),
      "1: foo, bar: NA, baz: ()"
    )
  })
})
