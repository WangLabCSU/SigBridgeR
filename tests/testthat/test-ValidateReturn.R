# ValidateReturn ----------------------------------------------------------

test_that("`ValidateReturn` accepts a list with `scRNA_data`", {
  expect_message(
    res <- ValidateReturn(compliant_screen),
    "Return value is a list with `scRNA_data` slot"
  )

  expect_null(res)
})

test_that("`ValidateReturn` detects a list without `scRNA_data`", {
  no_slot <- function(sc_data) {
    list(other = sc_data)
  }

  expect_message(
    res <- ValidateReturn(no_slot),
    "does not have `scRNA_data` slot"
  )

  expect_equal(res, list(error = 1L))
})

test_that("`ValidateReturn` detects a non-list return value", {
  not_list <- function(sc_data) {
    sc_data
  }

  expect_message(
    res <- ValidateReturn(not_list),
    "Return value is not a list"
  )

  expect_equal(res, list(error = 1L))
})

test_that("`ValidateReturn` handles an empty `list()`", {
  empty_list <- function(sc_data) {
    list()
  }

  expect_message(
    res <- ValidateReturn(empty_list),
    "does not have `scRNA_data` slot"
  )

  expect_equal(res, list(error = 1L))
})

test_that("`ValidateReturn` follows explicit `return()` calls", {
  explicit <- function(sc_data) {
    if (FALSE) {
      message("unreachable")
    }
    return(list(scRNA_data = sc_data))
  }

  expect_null(suppressMessages(ValidateReturn(explicit)))

  explicit_bad <- function(sc_data) {
    return(sc_data)
  }

  expect_equal(
    suppressMessages(ValidateReturn(explicit_bad)),
    list(error = 1L)
  )
})

test_that("`ValidateReturn` checks early `return()` calls", {
  early_good <- function(sc_data) {
    if (sc_data > 0) {
      return(list(scRNA_data = sc_data))
    }
    list(scRNA_data = sc_data)
  }

  early_bad <- function(sc_data) {
    if (sc_data > 0) {
      return(sc_data)
    }
    list(scRNA_data = sc_data)
  }

  expect_null(suppressMessages(ValidateReturn(early_good)))
  expect_equal(
    suppressMessages(ValidateReturn(early_bad)),
    list(error = 1L)
  )
})

test_that("`ValidateReturn` checks every `if` branch", {
  both_good <- function(sc_data, flag = TRUE) {
    if (flag) {
      list(scRNA_data = sc_data)
    } else {
      list(scRNA_data = sc_data)
    }
  }

  then_bad <- function(sc_data, flag = TRUE) {
    if (flag) {
      list(other = sc_data)
    } else {
      list(scRNA_data = sc_data)
    }
  }

  else_bad <- function(sc_data, flag = TRUE) {
    if (flag) {
      list(scRNA_data = sc_data)
    } else {
      list(other = sc_data)
    }
  }

  expect_null(suppressMessages(ValidateReturn(both_good)))

  expect_message(
    res <- ValidateReturn(then_bad),
    "does not have `scRNA_data` slot"
  )
  expect_equal(res, list(error = 1L))

  expect_message(
    res <- ValidateReturn(else_bad),
    "does not have `scRNA_data` slot"
  )
  expect_equal(res, list(error = 1L))
})

test_that("`ValidateReturn` treats a missing `else` as `NULL`", {
  no_else <- function(sc_data, flag = TRUE) {
    if (flag) {
      list(scRNA_data = sc_data)
    }
  }

  expect_message(
    res <- ValidateReturn(no_else),
    "Return value is not a list"
  )

  expect_equal(res, list(error = 1L))
})

# .vr_* helpers -----------------------------------------------------------

test_that("`.vr_collect_return_values` follows blocks, `if` and `return()`", {
  values <- .vr_collect_return_values(body(function(sc_data, flag = TRUE) {
    if (flag) {
      return(list(scRNA_data = sc_data))
    }
    list(scRNA_data = sc_data)
  }))

  expect_length(values, 2L)

  no_else <- .vr_collect_return_values(quote(if (flag) list(a = 1)))

  # the missing `else` is recorded as `NULL`
  expect_length(no_else, 2L)
  expect_null(no_else[[2L]])
})

test_that("`.vr_collect_return_values` keeps plain expressions as is", {
  expect_identical(
    .vr_collect_return_values(quote(sc_data)),
    list(quote(sc_data))
  )
  expect_identical(
    .vr_collect_return_values(quote(list(scRNA_data = 1L))),
    list(quote(list(scRNA_data = 1L)))
  )
})

test_that("`.vr_explicit_returns_in` skips nested functions", {
  nested <- .vr_collect_return_values(body(function(sc_data) {
    helper <- function() {
      return(1L)
    }
    list(scRNA_data = sc_data)
  }))

  expect_length(nested, 1L)
  expect_identical(nested[[1L]], quote(list(scRNA_data = sc_data)))
})

test_that("`.vr_is_list_creation` recognises `list()` calls only", {
  expect_true(.vr_is_list_creation(quote(list(a = 1))))
  expect_true(.vr_is_list_creation(quote(list())))
  expect_false(.vr_is_list_creation(quote(data.frame(a = 1))))
  expect_false(.vr_is_list_creation(quote(x)))
  expect_false(.vr_is_list_creation(NULL))
})

test_that("`.vr_list_has_scRNA_data` detects the slot", {
  expect_true(.vr_list_has_scRNA_data(quote(list(scRNA_data = 1))))
  expect_true(.vr_list_has_scRNA_data(quote(list(1, scRNA_data = 2))))
  expect_false(.vr_list_has_scRNA_data(quote(list(other = 1))))
  expect_false(.vr_list_has_scRNA_data(quote(list())))
  expect_false(.vr_list_has_scRNA_data(quote(scRNA_data)))
})

# safe_call_name / safe_call_args -----------------------------------------

test_that("`safe_call_name` tolerates non-call and formula input", {
  expect_identical(safe_call_name(quote(dir.create("x"))), "dir.create")
  expect_null(safe_call_name(quote(some_symbol)))
  expect_null(safe_call_name(~some_formula))
})

test_that("`safe_call_args` returns arguments or an empty list", {
  expect_length(safe_call_args(quote(f(1, 2))), 2L)
  expect_identical(safe_call_args(quote(some_symbol)), list())
  expect_identical(safe_call_args(~some_formula), list())
})
