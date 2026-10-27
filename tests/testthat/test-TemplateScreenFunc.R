# TemplateScreenFunc - argument validation -------------------------------

test_that("`TemplateScreenFunc` validates `filename`", {
  expect_error(TemplateScreenFunc(filename = 1), class = "chk_error")
  expect_error(TemplateScreenFunc(filename = character()), class = "chk_error")
  expect_error(
    TemplateScreenFunc(filename = c("a.R", "b.R")),
    class = "chk_error"
  )
  expect_error(
    TemplateScreenFunc(file.path(withr::local_tempdir(), "screen.txt")),
    "must have a '\\.R' extension"
  )
})

test_that("`TemplateScreenFunc` validates `func_name`", {
  path <- file.path(withr::local_tempdir(), "screen.R")

  expect_error(
    TemplateScreenFunc(path, func_name = "1bad", open = FALSE),
    "Invalid `func_name`"
  )
  expect_error(
    TemplateScreenFunc(path, func_name = "has space", open = FALSE),
    "Invalid `func_name`"
  )
  expect_error(
    TemplateScreenFunc(path, func_name = "a-b", open = FALSE),
    "Invalid `func_name`"
  )

  # nothing is written when validation fails
  expect_false(file.exists(path))
})

test_that("`TemplateScreenFunc` accepts valid `func_name` shapes", {
  path <- file.path(withr::local_tempdir(), "screen.R")

  expect_no_error(
    TemplateScreenFunc(path, func_name = "MyScreen_2.0", open = FALSE)
  )
  expect_true(any(grepl("^MyScreen_2.0 <- function\\(", readLines(path))))
})

# TemplateScreenFunc - writing --------------------------------------------

test_that("`TemplateScreenFunc` writes a template without documentation", {
  path <- file.path(withr::local_tempdir(), "my_screen.R")

  expect_message(
    res <- TemplateScreenFunc(path, func_name = "my_screen", open = FALSE),
    "successfully written to"
  )

  expect_identical(res, path)
  expect_true(file.exists(path))

  lines <- readLines(path)
  expect_false(any(grepl("^#'", lines)))
  expect_true(any(grepl("^my_screen <- function\\($", lines)))
  expect_false(any(grepl("FUNC_NAME", lines, fixed = TRUE)))
  expect_true(any(grepl("scRNA_data", lines)))
})

test_that("`TemplateScreenFunc` includes documentation when requested", {
  path <- file.path(withr::local_tempdir(), "documented.R")

  expect_message(
    TemplateScreenFunc(
      path,
      func_name = "documented_screen",
      documentation = TRUE,
      open = FALSE
    ),
    "Roxygen2 documentation included"
  )

  lines <- readLines(path)
  expect_true(any(grepl("^#' @title", lines)))
  expect_true(any(grepl("^#' @export", lines)))
  expect_true(any(grepl("^documented_screen <- function\\(", lines)))
})

test_that("`TemplateScreenFunc` returns the path invisibly", {
  path <- file.path(withr::local_tempdir(), "invisible.R")

  res <- withVisible(
    suppressMessages(TemplateScreenFunc(path, open = FALSE))
  )

  expect_false(res$visible)
  expect_identical(res$value, path)
})

test_that("`TemplateScreenFunc` uses `my_screen.R` by default", {
  withr::local_dir(withr::local_tempdir())

  expect_message(
    res <- TemplateScreenFunc(func_name = "default_screen", open = FALSE),
    "successfully written to"
  )

  expect_identical(res, "my_screen.R")
  expect_true(file.exists("my_screen.R"))
})

# TemplateScreenFunc - append / overwrite ---------------------------------

test_that("`TemplateScreenFunc` appends to an existing file by default", {
  path <- withr::local_tempfile(fileext = ".R", lines = "# my screens")

  suppressMessages(
    TemplateScreenFunc(path, func_name = "first_screen", open = FALSE)
  )
  suppressMessages(
    TemplateScreenFunc(path, func_name = "second_screen", open = FALSE)
  )

  lines <- readLines(path)
  expect_true("# my screens" %in% lines)
  expect_true(any(grepl("^first_screen <- function\\(", lines)))
  expect_true(any(grepl("^second_screen <- function\\(", lines)))

  # the appended template is separated from the previous content
  header_idx <- which(lines == "# Template screen function for SigBridgeR")
  expect_length(header_idx, 2L)
  expect_identical(lines[header_idx[2L] - 1L], "")
})

test_that("`TemplateScreenFunc` overwrites when the user agrees", {
  path <- withr::local_tempfile(fileext = ".R", lines = c("# old", "# old"))

  local_mocked_bindings(
    askYesNo = function(...) TRUE,
    .package = "utils"
  )

  expect_message(
    TemplateScreenFunc(
      path,
      func_name = "fresh_screen",
      append = FALSE,
      open = FALSE
    ),
    "successfully written to"
  )

  lines <- readLines(path)
  expect_false(any(grepl("^# old", lines)))
  expect_true(any(grepl("^fresh_screen <- function\\(", lines)))
})

test_that("`TemplateScreenFunc` aborts when the user declines overwrite", {
  path <- withr::local_tempfile(fileext = ".R", lines = "# keep me")

  local_mocked_bindings(
    askYesNo = function(...) FALSE,
    .package = "utils"
  )

  expect_error(
    TemplateScreenFunc(
      path,
      func_name = "nope_screen",
      append = FALSE,
      open = FALSE
    ),
    "Operation cancelled by user"
  )

  expect_identical(readLines(path), "# keep me")
})

test_that("`TemplateScreenFunc` skips confirmation for an empty file", {
  path <- withr::local_tempfile(fileext = ".R")

  local_mocked_bindings(
    askYesNo = function(...) stop("confirmation should not be requested"),
    .package = "utils"
  )

  expect_no_error(
    suppressMessages(
      TemplateScreenFunc(
        path,
        func_name = "empty_file",
        append = FALSE,
        open = FALSE
      )
    )
  )
  expect_true(any(grepl("^empty_file <- function\\(", readLines(path))))
})

# TemplateScreenFunc - `filename = "current"` -----------------------------

test_that("`TemplateScreenFunc` writes to the current script", {
  skip_if_not_installed("rstudioapi")

  current <- withr::local_tempfile(fileext = ".R", lines = "# current script")

  local_mocked_bindings(
    isAvailable = function() TRUE,
    getSourceEditorContext = function(...) list(path = current),
    .package = "rstudioapi"
  )

  expect_message(
    res <- TemplateScreenFunc(
      "current",
      func_name = "current_screen",
      open = FALSE
    ),
    "Using current script"
  )

  expect_identical(res, current)
  expect_true(any(grepl("^current_screen <- function\\(", readLines(current))))
})

test_that("`filename = 'current'` is case-insensitive", {
  skip_if_not_installed("rstudioapi")

  current <- withr::local_tempfile(fileext = ".R")

  local_mocked_bindings(
    isAvailable = function() TRUE,
    getSourceEditorContext = function(...) list(path = current),
    .package = "rstudioapi"
  )

  expect_no_error(
    suppressMessages(
      TemplateScreenFunc("CURRENT", func_name = "upper_screen", open = FALSE)
    )
  )
  expect_true(any(grepl("^upper_screen <- function\\(", readLines(current))))
})

test_that("`TemplateScreenFunc` aborts on `current` outside RStudio", {
  skip_if_not_installed("rstudioapi")

  local_mocked_bindings(
    isAvailable = function() FALSE,
    .package = "rstudioapi"
  )

  expect_error(
    TemplateScreenFunc("current", func_name = "no_session", open = FALSE),
    "requires RStudio or Positron"
  )
})

# TemplateScreenFunc - opening the file -----------------------------------

test_that("`TemplateScreenFunc` opens the file in the editor when available", {
  skip_if_not_installed("rstudioapi")

  path <- file.path(withr::local_tempdir(), "open_me.R")
  opened <- new.env(parent = emptyenv())

  local_mocked_bindings(
    isAvailable = function() TRUE,
    navigateToFile = function(file, ...) {
      opened$file <- file
      invisible(NULL)
    },
    .package = "rstudioapi"
  )

  expect_message(
    TemplateScreenFunc(path, func_name = "opened_screen", open = TRUE),
    "File opened in editor"
  )

  expect_identical(opened$file, path)
})

test_that("`TemplateScreenFunc` does not open the file when `open = FALSE`", {
  path <- file.path(withr::local_tempdir(), "closed.R")

  expect_no_error(
    suppressMessages(
      TemplateScreenFunc(path, func_name = "closed_screen", open = FALSE)
    )
  )
  expect_true(file.exists(path))
})

# build_template_body -----------------------------------------------------

test_that("`build_template_body` strips documentation on request", {
  body <- build_template_body("plain_screen", documentation = FALSE)

  expect_type(body, "character")
  expect_false(any(grepl("^#'", body)))
  expect_true(any(grepl("^plain_screen <- function\\(", body)))
  expect_false(any(grepl("FUNC_NAME", body, fixed = TRUE)))
})

test_that("`build_template_body` keeps documentation on request", {
  body <- build_template_body("documented_screen", documentation = TRUE)

  expect_true(any(grepl("^#' @title", body)))
  expect_true(any(grepl("^documented_screen <- function\\(", body)))
  expect_false(any(grepl("FUNC_NAME", body, fixed = TRUE)))
})

test_that("`build_template_body` keeps the expected return structure", {
  body <- trimws(build_template_body("screen_fun", documentation = FALSE))

  expect_in("scRNA_data = modified_sc_data", body)
  expect_in("sc_data,", body)
  expect_in("matched_bulk,", body)
  expect_in("phenotype,", body)
})
