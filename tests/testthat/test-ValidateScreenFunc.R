# ValidateScreenFunc - orchestration --------------------------------------

test_that("`ValidateScreenFunc` returns TRUE invisibly", {
  res <- withVisible(
    suppressMessages(ValidateScreenFunc(compliant_screen))
  )

  expect_false(res$visible)
  expect_true(res$value)
})

test_that("`ValidateScreenFunc` prints a header and a badge", {
  expect_message(
    ValidateScreenFunc(compliant_screen),
    "Screening Function Validation"
  )
  expect_message(ValidateScreenFunc(compliant_screen), "Start at")
  expect_message(ValidateScreenFunc(compliant_screen), "error")
})

test_that("`ValidateScreenFunc` reports problems without throwing", {
  incomplete <- function(x) {
    list(foo = x)
  }

  expect_no_error(
    res <- suppressMessages(ValidateScreenFunc(incomplete))
  )
  expect_true(res)

  expect_message(
    ValidateScreenFunc(incomplete),
    "Missing required arguments"
  )
})

test_that("`ValidateScreenFunc` skips the timer when tictoc is missing", {
  local_mocked_bindings(is_installed = function(...) FALSE)

  msgs <- suppressWarnings(
    capture_messages(ValidateScreenFunc(compliant_screen))
  )

  expect_false(any(grepl("Duration", msgs)))
})
