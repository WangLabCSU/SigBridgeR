# ValidateArgsInputs ------------------------------------------------------

test_that("`ValidateArgsInputs` accepts a fully specified interface", {
  expect_message(
    res <- ValidateArgsInputs(compliant_screen),
    "All input arguments explicitly specified"
  )

  expect_equal(res, list(error = 0L))
})

test_that("`ValidateArgsInputs` reports missing arguments", {
  incomplete <- function(sc_data, ...) {
    list(scRNA_data = sc_data)
  }

  expect_message(
    res <- ValidateArgsInputs(incomplete),
    "Missing required arguments"
  )

  expect_equal(res, list(error = 1L))
})

test_that("`ValidateArgsInputs` explains every missing argument", {
  no_args <- function() NULL

  expect_message(
    ValidateArgsInputs(no_args),
    "More arguments should be added"
  )
  expect_message(ValidateArgsInputs(no_args), "phenotype_class")
})

test_that("`ValidateArgsInputs` accepts a partial interface", {
  partial <- function(sc_data, matched_bulk, phenotype) {
    list(scRNA_data = sc_data)
  }

  expect_message(
    res <- ValidateArgsInputs(partial),
    "Missing required arguments"
  )
  expect_message(ValidateArgsInputs(partial), "label_type")

  expect_equal(res, list(error = 1L))
})
