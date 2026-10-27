# ValidateDirName ---------------------------------------------------------

test_that("`ValidateDirName` is silent without `dir.create()`", {
  expect_null(suppressMessages(ValidateDirName(compliant_screen)))
})

test_that("`ValidateDirName` accepts folders ending with `_res`", {
  with_suffix <- function(sc_data) {
    dir.create("custom_res")
    list(scRNA_data = sc_data)
  }

  expect_message(
    res <- ValidateDirName(with_suffix),
    "ends with '_res'"
  )

  expect_null(res)
})

test_that("`ValidateDirName` notes folders without `_res`", {
  no_suffix <- function(sc_data) {
    dir.create("results")
    list(scRNA_data = sc_data)
  }

  expect_message(
    res <- ValidateDirName(no_suffix),
    "Folder name does not end with '_res'"
  )
  expect_message(ValidateDirName(no_suffix), "'results'")
  # the suggestion is derived from the current folder name
  expect_message(ValidateDirName(no_suffix), "'results_res'")

  expect_equal(res, list(note = 1L))
})

test_that("`ValidateDirName` matches the suffix at the end only", {
  contains_suffix <- function(sc_data) {
    dir.create("my_results")
    list(scRNA_data = sc_data)
  }

  expect_message(
    res <- ValidateDirName(contains_suffix),
    "Folder name does not end with '_res'"
  )

  expect_equal(res, list(note = 1L))
})

test_that("`ValidateDirName` ignores trailing separators", {
  trailing_slash <- function(sc_data) {
    dir.create("slashed_res/")
    list(scRNA_data = sc_data)
  }

  expect_null(suppressMessages(ValidateDirName(trailing_slash)))
})

test_that("`ValidateDirName` locates the `path` argument", {
  named_path <- function(sc_data) {
    dir.create(showWarnings = FALSE, path = "named_res")
    list(scRNA_data = sc_data)
  }

  positional_path <- function(sc_data) {
    dir.create(showWarnings = FALSE, "positional_res")
    list(scRNA_data = sc_data)
  }

  expect_null(suppressMessages(ValidateDirName(named_path)))
  expect_null(suppressMessages(ValidateDirName(positional_path)))
})

test_that("`ValidateDirName` warns on `dir.create()` without `path`", {
  no_path <- function(sc_data) {
    dir.create()
    list(scRNA_data = sc_data)
  }

  expect_warning(
    res <- ValidateDirName(no_path),
    "Failed to get folder name"
  )

  expect_null(res)
})

test_that("`ValidateDirName` resolves folder names held in variables", {
  good_var <- function(sc_data) {
    out <- "variable_res"
    dir.create(out)
    list(scRNA_data = sc_data)
  }

  bad_var <- function(sc_data) {
    out <- "variable"
    dir.create(out)
    list(scRNA_data = sc_data)
  }

  expect_null(suppressMessages(ValidateDirName(good_var)))
  expect_equal(
    suppressMessages(ValidateDirName(bad_var)),
    list(note = 1L)
  )
})

test_that("`ValidateDirName` counts every non-compliant folder", {
  mixed <- function(sc_data) {
    dir.create("good_res")
    dir.create("bad_one")
    dir.create("bad_two")
    list(scRNA_data = sc_data)
  }

  expect_equal(
    suppressMessages(ValidateDirName(mixed)),
    list(note = 2L)
  )
})

test_that("`ValidateDirName` warns on unevaluable folder names", {
  complex_path <- function(sc_data) {
    dir.create(file.path("a", "b_res"))
    list(scRNA_data = sc_data)
  }

  expect_warning(
    res <- ValidateDirName(complex_path),
    "Failed to get folder name"
  )

  expect_null(res)
})

test_that("`ValidateDirName` warns on unknown folder variables", {
  unknown_var <- function(sc_data) {
    dir.create(undefined_folder_var)
    list(scRNA_data = sc_data)
  }

  # `get_var_value()` aborts on symbols it cannot resolve
  expect_warning(
    res <- ValidateDirName(unknown_var),
    "Failed to get folder name"
  )

  expect_null(res)
})

test_that("`ValidateDirName` reports a folder created twice only once", {
  duplicated <- function(sc_data, flag = TRUE) {
    if (flag) {
      dir.create("same_folder")
    } else {
      dir.create("same_folder")
    }
    list(scRNA_data = sc_data)
  }

  expect_equal(
    suppressMessages(ValidateDirName(duplicated)),
    list(note = 1L)
  )
})

test_that("`ValidateDirName` skips expressions that are not inspectable", {
  with_formula <- function(sc_data) {
    ~sc_data
    dir.create("formula_res")
    list(scRNA_data = sc_data)
  }

  # `is_call()` accepts formulas while `call_name()` rejects them
  expect_no_error(suppressMessages(ValidateDirName(with_formula)))
  expect_null(suppressMessages(ValidateDirName(with_formula)))
})

# .vdn_* helpers ----------------------------------------------------------

test_that("`.vdn_find_dir_create_calls` collects nested calls", {
  body <- body(function(sc_data) {
    for (i in 1:2) {
      dir.create(paste0("folder", i))
    }
    dir.create("second_res")
  })

  calls <- .vdn_find_dir_create_calls(body)

  expect_length(calls, 2L)
  expect_identical(safe_call_name(calls[[1L]]), "dir.create")
})

test_that("`.vdn_find_dir_create_calls` tolerates formulas and literals", {
  expect_identical(.vdn_find_dir_create_calls(quote(~x)), list())
  expect_identical(.vdn_find_dir_create_calls(1L), list())
})

test_that("`.vdn_path_arg` matches `path` by name or position", {
  expect_identical(
    .vdn_path_arg(quote(dir.create("plain_res"))),
    "plain_res"
  )
  expect_identical(
    .vdn_path_arg(quote(dir.create(recursive = TRUE, "positional_res"))),
    "positional_res"
  )
  expect_identical(
    .vdn_path_arg(quote(dir.create(recursive = TRUE, path = "named_res"))),
    "named_res"
  )
  expect_null(.vdn_path_arg(quote(dir.create())))
  expect_null(.vdn_path_arg(quote(dir.create(recursive = TRUE))))
})

test_that("`.vdn_folder_name` resolves literals, variables and failures", {
  runner <- function(sc_data) {
    out <- "variable_res"
    dir.create(out)
    list(scRNA_data = sc_data)
  }

  expect_identical(.vdn_folder_name("literal_res", runner), "literal_res")
  expect_identical(.vdn_folder_name(quote(out), runner), "variable_res")
  expect_identical(.vdn_folder_name(123, runner), "123")
  expect_identical(.vdn_folder_name("slashed_res/", runner), "slashed_res")

  # unevaluable inputs are reported as `NULL`
  expect_null(.vdn_folder_name(quote(missing_var), runner))
  expect_null(.vdn_folder_name(quote(file.path("a", "b")), runner))
  expect_null(.vdn_folder_name(NA_character_, runner))
  expect_null(.vdn_folder_name(character(), runner))
  expect_null(.vdn_folder_name(quote(x), function(x) list(scRNA_data = x)))

  vector_var <- function(sc_data) {
    out <- c("a_res", "b_res")
    dir.create(out)
    list(scRNA_data = sc_data)
  }
  expect_null(.vdn_folder_name(quote(out), vector_var))
})

test_that("`.vdn_has_res_suffix` only matches the end of the name", {
  expect_true(.vdn_has_res_suffix("custom_res"))
  expect_true(.vdn_has_res_suffix("_res"))
  expect_false(.vdn_has_res_suffix("my_results"))
  expect_false(.vdn_has_res_suffix("res"))
  expect_false(.vdn_has_res_suffix("results_res_backup"))
})
