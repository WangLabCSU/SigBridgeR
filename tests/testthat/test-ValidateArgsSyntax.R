# ValidateArgsSyntax ------------------------------------------------------

test_that("`ValidateArgsSyntax` passes a clean function", {
  clean <- function(sc_data, matched_bulk) {
    list(scRNA_data = sc_data, bulk = matched_bulk)
  }

  skip_if_not_installed("codetools")

  expect_message(
    res <- ValidateArgsSyntax(clean),
    "Syntax check passed"
  )

  expect_equal(res, list(error = 0L, warn = 0L))
})

test_that("`ValidateArgsSyntax` flags unused parameters", {
  unused <- function(sc_data, matched_bulk) {
    list(scRNA_data = sc_data)
  }

  skip_if_not_installed("codetools")

  expect_message(
    res <- ValidateArgsSyntax(unused),
    "Syntax error in function"
  )

  expect_equal(res$error, 1L)
})

test_that("`ValidateArgsSyntax` flags suspicious globals", {
  bad_global <- function(sc_data) {
    list(scRNA_data = sc_data, extra = undefined_object_xyz)
  }

  skip_if_not_installed("codetools")

  expect_message(
    res <- ValidateArgsSyntax(bad_global),
    "Suspicious global variables"
  )

  expect_equal(res$warn, 1L)
})

test_that("`ValidateArgsSyntax` warns when codetools is unavailable", {
  local_mocked_bindings(is_installed = function(...) FALSE)

  expect_warning(
    ValidateArgsSyntax(compliant_screen),
    "not supported"
  )
})
