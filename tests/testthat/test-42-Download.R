# Tests for fileDownload() (R/42-Download.R) --------------------------------

test_that("fileDownload creates the destination directory and writes the file", {
  root <- withr::local_tempdir()
  dest <- file.path(root, "nested", "dir", "data.rds")

  local_mocked_bindings(
    curl_download = function(url, destfile, ...) {
      writeLines("content", destfile)
      TRUE
    },
    .package = "curl"
  )

  expect_invisible(
    res <- fileDownload("https://example.test/data.rds", dest, quiet = TRUE, retries = 1L)
  )
  expect_identical(res, dest)
  expect_true(dir.exists(dirname(dest)))
  expect_identical(readLines(dest), "content")
  expect_false(file.exists(paste0(dest, ".tmp")))
})

test_that("fileDownload passes timeout and handle settings to curl", {
  root <- withr::local_tempdir()
  dest <- file.path(root, "data.rds")
  handle_args <- NULL
  download_args <- NULL

  local_mocked_bindings(
    new_handle = function(...) {
      handle_args <<- list(...)
      NULL
    },
    curl_download = function(url, destfile, quiet, mode, handle, ...) {
      download_args <<- list(
        url = url,
        destfile = destfile,
        quiet = quiet,
        mode = mode,
        handle = handle
      )
      writeLines("x", destfile)
      TRUE
    },
    .package = "curl"
  )

  fileDownload("https://example.test/f", dest, quiet = TRUE, retries = 1L, timeout = 42L)

  expect_identical(handle_args$followlocation, TRUE)
  expect_identical(handle_args$connecttimeout, 20L)
  expect_identical(handle_args$timeout, 42L)
  expect_match(handle_args$useragent, "^R/")
  expect_identical(download_args$url, "https://example.test/f")
  expect_identical(download_args$destfile, paste0(dest, ".tmp"))
  expect_true(download_args$quiet)
  expect_identical(download_args$mode, "wb")
  expect_null(download_args$handle)
})

test_that("fileDownload retries a failed download", {
  root <- withr::local_tempdir()
  dest <- file.path(root, "data.rds")
  n <- 0L

  local_mocked_bindings(
    curl_download = function(url, destfile, ...) {
      n <<- n + 1L
      if (n < 3L) {
        stop("temporary failure")
      }
      writeLines("ok", destfile)
      TRUE
    },
    .package = "curl"
  )
  local_mocked_bindings(Sys.sleep = function(...) NULL, .package = "base")

  fileDownload("https://example.test/f", dest, quiet = TRUE, retries = 3L)

  expect_identical(n, 3L)
  expect_identical(readLines(dest), "ok")
  expect_false(file.exists(paste0(dest, ".tmp")))
})

test_that("fileDownload falls back to copy when rename fails", {
  root <- withr::local_tempdir()
  dest <- file.path(root, "data.rds")

  local_mocked_bindings(
    curl_download = function(url, destfile, ...) {
      writeLines("copied", destfile)
      TRUE
    },
    .package = "curl"
  )
  local_mocked_bindings(file.rename = function(from, to) FALSE, .package = "base")

  fileDownload("https://example.test/f", dest, quiet = TRUE, retries = 1L)

  expect_identical(readLines(dest), "copied")
  expect_false(file.exists(paste0(dest, ".tmp")))
})

test_that("fileDownload stops after all retries fail", {
  root <- withr::local_tempdir()
  dest <- file.path(root, "data.rds")
  n <- 0L

  local_mocked_bindings(
    curl_download = function(...) {
      n <<- n + 1L
      stop("connection refused")
    },
    .package = "curl"
  )
  local_mocked_bindings(Sys.sleep = function(...) NULL, .package = "base")

  expect_error(
    fileDownload("https://example.test/f", dest, quiet = TRUE, retries = 3L),
    "Download failed after 3 attempts: connection refused"
  )
  expect_identical(n, 3L)
  expect_false(file.exists(dest))
  expect_false(file.exists(paste0(dest, ".tmp")))
})
