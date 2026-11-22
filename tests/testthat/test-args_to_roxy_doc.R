# args_to_roxy_doc ---------------------------------------------------------------

a2r_parse_rd <- function(x) {
  con <- textConnection(x)
  on.exit(close(con))
  tools::parse_Rd(con)
}

a2r_make_dev_pkg <- function(dir) {
  dir.create(file.path(dir, "man"))
  writeLines(
    c("Package: devpkg", "Version: 0.0.1"),
    file.path(dir, "DESCRIPTION")
  )
  writeLines(
    c(
      "\\name{foo}",
      "\\alias{foo}",
      "\\arguments{",
      "\\item{x}{X description.}",
      "\\item{y, z}{Shared description.}",
      "}"
    ),
    file.path(dir, "man", "foo.Rd")
  )
  invisible(dir)
}

test_that("args_to_roxy_doc reads arguments from a development package", {
  dir <- withr::local_tempdir()
  a2r_make_dev_pkg(dir)

  expect_equal(
    args_to_roxy_doc("foo", path = dir),
    c(
      "- `x`: X description.",
      "- `y`: Shared description.",
      "- `z`: Shared description."
    )
  )
})

test_that("args_to_roxy_doc matches an Rd topic by \\name", {
  dir <- withr::local_tempdir()
  a2r_make_dev_pkg(dir)
  writeLines(
    c("\\name{only_name}", "\\arguments{\\item{a}{A description.}}"),
    file.path(dir, "man", "only_name.Rd")
  )

  expect_equal(args_to_roxy_doc("only_name", path = dir), "- `a`: A description.")
})

test_that("args_to_roxy_doc returns character(0) when a topic has no arguments", {
  dir <- withr::local_tempdir()
  a2r_make_dev_pkg(dir)
  writeLines(
    c("\\name{empty}", "\\alias{empty}", "\\description{No arguments here.}"),
    file.path(dir, "man", "empty.Rd")
  )

  expect_identical(args_to_roxy_doc("empty", path = dir), character(0L))
})

test_that("args_to_roxy_doc reads documentation from an installed package", {
  res <- args_to_roxy_doc("lm", package = "stats")

  expect_match(res, "^- `", all = TRUE)
  expect_match(res[[1L]], "formula")
})

test_that("args_to_roxy_doc infers the package from a function or a qualified name", {
  from_function <- args_to_roxy_doc(stats::sd)
  from_string <- args_to_roxy_doc("stats::sd")

  expect_identical(from_function, from_string)
  expect_match(from_function[[1L]], "^- `x`:")
})

test_that("args_to_roxy_doc rejects `fn` values that are not a single function or string", {
  expect_error(args_to_roxy_doc(1L), "must be a function or a single string")
  expect_error(args_to_roxy_doc(c("a", "b")), "must be a function or a single string")
})

test_that("args_to_roxy_doc warns when the package cannot be determined", {
  expect_warning(
    res <- args_to_roxy_doc("zzz_no_such_fn_xyz"),
    "Cannot determine the package"
  )
  expect_identical(res, character(0L))
})

test_that("args_to_roxy_doc warns when the requested package is not installed", {
  expect_warning(
    res <- args_to_roxy_doc("foo", package = "NotARealPkgXYZ"),
    "is not installed"
  )
  expect_identical(res, character(0L))
})

test_that("args_to_roxy_doc warns when no Rd topic matches", {
  expect_warning(
    res <- args_to_roxy_doc("zzz_no_such_topic", package = "stats"),
    "No Rd documentation found"
  )
  expect_identical(res, character(0L))
})

# find_rd_by_alias --------------------------------------------------------

test_that("find_rd_by_alias matches \\alias and \\name tags", {
  rd <- a2r_parse_rd("\\name{longName}\\alias{shortName}\\arguments{}")

  expect_identical(find_rd_by_alias(list(rd), "shortName"), rd)
  expect_identical(find_rd_by_alias(list(rd), "longName"), rd)
  expect_null(find_rd_by_alias(list(rd), "missing"))
})

# rd_arguments_to_md ------------------------------------------------------

test_that("rd_arguments_to_md collapses whitespace and splits multi-argument items", {
  rd <- a2r_parse_rd("\\arguments{\\item{x, y}{First   line.\n Second line.}}")

  expect_equal(
    rd_arguments_to_md(rd),
    c("- `x`: First line. Second line.", "- `y`: First line. Second line.")
  )
})

test_that("rd_arguments_to_md returns character(0) without an arguments section", {
  expect_identical(
    rd_arguments_to_md(a2r_parse_rd("\\name{foo}\\description{None.}")),
    character(0L)
  )
})

# rd_to_md ----------------------------------------------------------------

test_that("rd_to_md converts inline Rd markup to Markdown", {
  rd <- a2r_parse_rd(paste0(
    "\\arguments{\\item{x}{",
    "\\verb{v} \\samp{s} \\env{e} \\file{f} \\pkg{p} \\command{c} \\option{o} ",
    "\\emph{e2} \\strong{b2} \\ldots}}"
  ))

  expect_equal(
    rd_arguments_to_md(rd),
    "- `x`: `v` `s` `e` `f` `p` `c` `o` *e2* **b2** ..."
  )
})

test_that("rd_to_md converts links, remaining tags, and plain text", {
  rd <- a2r_parse_rd("\\arguments{\\item{x}{\\href{https://ex.com}{link}}}")
  expect_equal(rd_arguments_to_md(rd), "- `x`: [link](https://ex.com)")

  expect_equal(rd_to_md(structure(list("txt"), Rd_tag = "\\em")), "*txt*")
  expect_equal(rd_to_md(structure(list("txt"), Rd_tag = "\\bold")), "**txt**")
  expect_equal(rd_to_md(structure(list("txt"), Rd_tag = "\\item")), " txt")
  expect_equal(rd_to_md("plain"), "plain")
  expect_equal(rd_to_md(structure("...", Rd_tag = "\\dots")), "...")
})
