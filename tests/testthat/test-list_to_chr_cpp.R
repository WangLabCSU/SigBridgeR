# Tests for src/list_to_chr.cpp — list_to_chr() and list_to_chr_formated(),
# exercised through their R (Rcpp-exported) interfaces and the public wrappers.
#
# Run independently with:
#   devtools::load_all()
#   testthat::test_file("tests/testthat/test-list_to_chr_cpp.R")

describe("list_to_chr - list rendering", {
  it("renders a flat named list", {
    expect_equal(
      list_to_chr(list(a = 1, b = 2)),
      "a: 1, b: 2"
    )
  })

  it("uses 1-based indices for unnamed elements", {
    expect_equal(
      list_to_chr(list(10, 20)),
      "1: 10, 2: 20"
    )
  })

  it("renders nested lists with parentheses", {
    expect_equal(
      list_to_chr(list(a = list(b = 1, c = 2))),
      "a: (b: 1, c: 2)"
    )
  })

  it("handles deeply nested lists", {
    expect_equal(
      list_to_chr(list(a = list(b = list(c = 3)))),
      "a: (b: (c: 3))"
    )
  })

  it("mixes named and unnamed elements", {
    expect_equal(
      list_to_chr(list(a = 1, 2, c = 3)),
      "a: 1, 2: 2, c: 3"
    )
  })

  it("returns an empty string for an empty list", {
    expect_equal(list_to_chr(list()), "")
  })

  it("uses index for empty-string names", {
    x <- setNames(list(1, 2), c("", "b"))
    expect_equal(list_to_chr(x), "1: 1, b: 2")
  })

  it("uses index for NA names", {
    x <- setNames(list(1, 2), c(NA_character_, "b"))
    expect_equal(list_to_chr(x), "1: 1, b: 2")
  })
})

describe("list_to_chr - value types", {
  it("renders NA as NA", {
    expect_equal(list_to_chr(list(x = NA)), "x: NA")
  })

  it("renders NULL as () (Rf_isNewList(NULL) is TRUE on modern R)", {
    expect_equal(list_to_chr(list(x = NULL)), "x: ()")
  })

  it("renders atomic vectors as comma-separated values", {
    expect_equal(
      list_to_chr(list(x = c(1, 2, 3))),
      "x: 1, 2, 3"
    )
  })

  it("renders character vectors containing NA", {
    expect_equal(
      list_to_chr(list(x = c("a", NA, "c"))),
      "x: a, NA, c"
    )
  })

  it("renders logical vectors", {
    expect_equal(
      list_to_chr(list(x = c(TRUE, FALSE, NA))),
      "x: TRUE, FALSE, NA"
    )
  })

  it("renders functions as <function>", {
    expect_equal(
      list_to_chr(list(f = function(x) x)),
      "f: <function>"
    )
  })

  it("renders environments as <environment>", {
    expect_equal(
      list_to_chr(list(e = new.env())),
      "e: <environment>"
    )
  })
})

describe("list_to_chr - public interface", {
  it("returns a single character string matching the C++ wrapper", {
    x <- list(a = 1, b = list(c = 2))
    expect_identical(list_to_chr(x), list_to_chr(x))
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

# ---- list_to_chr_formated(): one indented line per element ------------------

describe("list_to_chr_formated - list rendering", {
  it("renders a flat named list, one element per line", {
    expect_equal(
      list_to_chr_formated(list(a = 1, b = 2)),
      "a: 1\nb: 2"
    )
  })

  it("uses 1-based indices for unnamed elements", {
    expect_equal(
      list_to_chr_formated(list(10, 20)),
      "1: 10\n2: 20"
    )
  })

  it("mixes named and unnamed elements", {
    expect_equal(
      list_to_chr_formated(list(a = 1, 2, c = 3)),
      "a: 1\n2: 2\nc: 3"
    )
  })

  it("uses index for empty-string names", {
    x <- setNames(list(1, 2), c("", "b"))
    expect_equal(list_to_chr_formated(x), "1: 1\nb: 2")
  })

  it("uses index for NA names", {
    x <- setNames(list(1, 2), c(NA_character_, "b"))
    expect_equal(list_to_chr_formated(x), "1: 1\nb: 2")
  })

  it("returns an empty string for an empty list", {
    expect_equal(list_to_chr_formated(list()), "")
  })

  it("has no trailing newline", {
    out <- list_to_chr_formated(list(a = 1, b = 2))
    expect_false(grepl("\n$", out))
  })

  it("emits exactly one line per top-level element", {
    out <- list_to_chr_formated(list(a = 1, b = 2, c = 3))
    expect_length(strsplit(out, "\n", fixed = TRUE)[[1]], 3L)
  })
})

describe("list_to_chr_formated - nesting and indentation", {
  it("puts a nested list on its own indented lines", {
    expect_equal(
      list_to_chr_formated(list(a = 1, b = list(c = 2, d = 3))),
      "a: 1\nb:\n  c: 2\n  d: 3"
    )
  })

  it("indents each deeper level by two more spaces", {
    expect_equal(
      list_to_chr_formated(list(a = list(b = list(c = 3)))),
      "a:\n  b:\n    c: 3"
    )
  })

  it("keeps the documented multi-line format", {
    x <- list(
      a = 1,
      b = list(c = 2, d = list(e = "hello")),
      f = c(1, 2, 3)
    )
    expect_equal(
      list_to_chr_formated(x),
      "a: 1\nb:\n  c: 2\n  d:\n    e: hello\nf: 1, 2, 3"
    )
  })

  it("resumes the parent indentation after a nested block", {
    expect_equal(
      list_to_chr_formated(list(a = list(b = 1), c = 2)),
      "a:\n  b: 1\nc: 2"
    )
  })

  it("uses positional keys for unnamed nested elements", {
    expect_equal(
      list_to_chr_formated(list(a = list(1, 2))),
      "a:\n  1: 1\n  2: 2"
    )
  })

  it("renders an empty nested list inline as ()", {
    expect_equal(
      list_to_chr_formated(list(a = 1, x = list(), b = 2)),
      "a: 1\nx: ()\nb: 2"
    )
  })

  it("expands a classed list such as a data frame as a nested list", {
    x <- list(df = data.frame(x = 1:2, y = c("a", "b")))
    expect_equal(
      list_to_chr_formated(x),
      "df:\n  x: 1, 2\n  y: a, b"
    )
  })
})

describe("list_to_chr_formated - value types", {
  it("renders NA as NA", {
    expect_equal(list_to_chr_formated(list(x = NA)), "x: NA")
  })

  it("renders NULL as NULL", {
    expect_equal(list_to_chr_formated(list(x = NULL)), "x: NULL")
  })

  it("renders atomic vectors on a single comma-separated line", {
    expect_equal(
      list_to_chr_formated(list(x = c(1, 2, 3))),
      "x: 1, 2, 3"
    )
  })

  it("renders character vectors containing NA", {
    expect_equal(
      list_to_chr_formated(list(x = c("a", NA, "c"))),
      "x: a, NA, c"
    )
  })

  it("renders logical vectors", {
    expect_equal(
      list_to_chr_formated(list(x = c(TRUE, FALSE, NA))),
      "x: TRUE, FALSE, NA"
    )
  })

  it("renders functions as <function>", {
    expect_equal(
      list_to_chr_formated(list(f = function(x) x)),
      "f: <function>"
    )
  })

  it("renders environments as <environment>", {
    expect_equal(
      list_to_chr_formated(list(e = new.env())),
      "e: <environment>"
    )
  })
})

describe("list_to_chr_formated - public interface", {
  it("returns a single character string", {
    x <- list(a = 1, b = list(c = 2))
    expect_type(list_to_chr_formated(x), "character")
    expect_length(list_to_chr_formated(x), 1)
  })

  it("is idempotent", {
    x <- list(a = 1, b = list(c = 2))
    expect_identical(list_to_chr_formated(x), list_to_chr_formated(x))
  })

  it("matches list_to_chr() for a flat list when lines are joined", {
    x <- list(a = 1, b = 2, c = c(3, 4))
    expect_identical(
      gsub("\n", ", ", list_to_chr_formated(x)),
      list_to_chr(x)
    )
  })
})
