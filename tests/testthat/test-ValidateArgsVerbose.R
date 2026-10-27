# ValidateArgsVerbose -----------------------------------------------------

test_that("`ValidateArgsVerbose` notes when verbose is not supported", {
  no_verbose <- function(sc_data) {
    list(scRNA_data = sc_data)
  }

  expect_message(
    res <- ValidateArgsVerbose(no_verbose),
    "Verbose control not supported"
  )

  expect_equal(res, list(note = 1L))
})

test_that("`ValidateArgsVerbose` warns when `verbose` is ignored", {
  dots_only <- function(sc_data, ...) {
    list(scRNA_data = sc_data)
  }

  expect_message(
    res <- ValidateArgsVerbose(dots_only),
    "Verbose control supported but not used"
  )

  expect_equal(res, list(warn = 1L))
})

test_that("`ValidateArgsVerbose` warns when a `verbose` argument is unused", {
  unused_verbose <- function(sc_data, verbose) {
    list(scRNA_data = sc_data)
  }

  expect_message(
    res <- ValidateArgsVerbose(unused_verbose),
    "Verbose control supported but not used"
  )

  expect_equal(res, list(warn = 1L))
})

test_that("`ValidateArgsVerbose` accepts a used `verbose` argument", {
  expect_message(
    res <- ValidateArgsVerbose(compliant_screen),
    "Verbose control supported"
  )

  expect_null(res)
})
