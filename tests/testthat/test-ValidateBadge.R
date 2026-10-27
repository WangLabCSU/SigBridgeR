# ValidateBadge -----------------------------------------------------------

test_that("`ValidateBadge` merges and pads check results", {
  expect_message(
    ValidateBadge(list(error = 1L), list(warn = 2L), list(note = 1L)),
    "1 error"
  )
  expect_message(
    ValidateBadge(list(error = 1L), list(error = 2L)),
    "3 errors"
  )
  expect_message(ValidateBadge(list(error = 0L)), "0 error")
})

test_that("`ValidateBadge` pluralises labels", {
  expect_message(ValidateBadge(list(warn = 1L)), "1 warning")
  expect_message(ValidateBadge(list(warn = 2L)), "2 warnings")
  expect_message(ValidateBadge(list(note = 1L)), "1 note")
  expect_message(ValidateBadge(list(note = 2L)), "2 notes")
})
