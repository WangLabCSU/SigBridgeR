# validate_* helpers ------------------------------------------------------

test_that("`validate_*` helpers tag their messages", {
  expect_message(validate_error("boom"), "\\.\\.\\. ERROR")
  expect_message(validate_warn("boom"), "\\.\\.\\. WARNING")
  expect_message(validate_note("boom"), "\\.\\.\\. NOTE")
  expect_message(validate_success("all good"), "all good")
})

test_that("`validate_explain` reports every message", {
  expect_message(validate_explain("first"), "first")
  expect_message(validate_explain(c("first", "second")), "first")
  expect_message(validate_explain(c("first", "second")), "second")
})

test_that("`validate_explain_speaker` formats single line messages", {
  expect_message(
    validate_explain_speaker("single line", "scRNA_data", "error"),
    "scRNA_data: single line"
  )
})

test_that("`validate_explain_speaker` lists additional lines", {
  messages <- c("head line", "bullet one", "bullet two")

  expect_message(
    validate_explain_speaker(messages, "scRNA_data", "warn"),
    "scRNA_data: head line"
  )
  expect_message(
    validate_explain_speaker(messages, "scRNA_data", "note"),
    "bullet two"
  )
})

test_that("`validate_explain_speaker` rejects unknown types", {
  expect_error(
    validate_explain_speaker("message", "scRNA_data", "invalid"),
    "should be one of"
  )
})
