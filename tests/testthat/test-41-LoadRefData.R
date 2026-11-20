# Tests for LoadRefData() (R/41-LoadRefData.R) -------------------------------

write_ref_cache <- function(path, data_type, value) {
  saveRDS(value, file.path(path, paste0(data_type, "_ref_data.rds")))
  invisible(path)
}

test_that("LoadRefData validates timeout", {
  tmp <- withr::local_tempdir()

  expect_error(LoadRefData("survival", path = tmp, timeout = 1.5), class = "chk_error")
  expect_error(LoadRefData("survival", path = tmp, timeout = "60"), class = "chk_error")
})

test_that("LoadRefData rejects an unknown data type", {
  tmp <- withr::local_tempdir()
  expect_error(LoadRefData("bogus", path = tmp), class = "MatchArgError")
})

test_that("LoadRefData requires path to be an existing directory", {
  expect_error(
    LoadRefData("survival", path = "/this/path/does/not/exist"),
    class = "chk_error"
  )

  file_path <- withr::local_tempfile()
  expect_error(LoadRefData("survival", path = file_path), class = "chk_error")
})

test_that("LoadRefData reads and returns cached data", {
  tmp <- withr::local_tempdir()
  value <- list(source = "cache", n = 3L)
  write_ref_cache(tmp, "survival", value)

  expect_message(
    res <- LoadRefData("survival", path = tmp, verbose = TRUE),
    "Found cached data"
  )
  expect_identical(res, value)
  expect_no_message(LoadRefData("survival", path = tmp, verbose = FALSE))
})

test_that("LoadRefData resolves data_type by exact and partial match", {
  tmp <- withr::local_tempdir()
  for (dt in c("survival", "binary", "continuous")) {
    write_ref_cache(tmp, dt, dt)
  }

  expect_identical(LoadRefData("survival", path = tmp, verbose = FALSE), "survival")
  expect_identical(LoadRefData("bin", path = tmp, verbose = FALSE), "binary")
  expect_identical(LoadRefData("cont", path = tmp, verbose = FALSE), "continuous")
})

test_that("LoadRefData downloads and caches missing data", {
  tmp <- withr::local_tempdir()
  value <- list(source = "download")
  calls <- list()

  local_mocked_bindings(
    fileDownload = function(url, destfile, mode, quiet, ...) {
      calls[[length(calls) + 1L]] <<- list(
        url = url,
        destfile = destfile,
        mode = mode,
        quiet = quiet
      )
      saveRDS(value, destfile)
      invisible(destfile)
    }
  )

  res <- LoadRefData("survival", path = tmp, timeout = 60L, verbose = FALSE)

  expect_identical(res, value)
  expect_length(calls, 1L)
  expect_match(calls[[1L]]$url, "survival_example_data\\.rds$")
  expect_identical(calls[[1L]]$mode, "wb")
  expect_false(calls[[1L]]$quiet)
  expect_identical(calls[[1L]]$destfile, file.path(tmp, "survival_ref_data.rds"))
  expect_true(file.exists(calls[[1L]]$destfile))
})

test_that("LoadRefData falls back to the next download source", {
  tmp <- withr::local_tempdir()
  urls <- character()

  local_mocked_bindings(
    fileDownload = function(url, destfile, ...) {
      urls <<- c(urls, url)
      if (length(urls) == 1L) {
        writeLines("partial", destfile)
        stop("first source unavailable")
      }
      saveRDS(list(source = "fallback"), destfile)
      invisible(destfile)
    }
  )

  res <- LoadRefData("binary", path = tmp, verbose = FALSE)

  expect_identical(res, list(source = "fallback"))
  expect_length(urls, 2L)
  expect_match(urls[[1L]], "raw.githubusercontent.com")
  expect_match(urls[[2L]], "raw.githubusercontent.com")
  expect_false(identical(urls[[1L]], urls[[2L]]))
  expect_match(urls[[2L]], "binary_example_data\\.rds$")
})

test_that("LoadRefData aborts after all download sources fail", {
  tmp <- withr::local_tempdir()
  n <- 0L

  local_mocked_bindings(
    fileDownload = function(url, destfile, ...) {
      n <<- n + 1L
      writeLines("partial", destfile)
      stop("network down")
    }
  )

  expect_error(
    LoadRefData("continuous", path = tmp, verbose = FALSE),
    "All download attempts failed"
  )
  expect_identical(n, 3L)
  expect_false(file.exists(file.path(tmp, "continuous_ref_data.rds")))
})

test_that("LoadRefData reports progress when verbose is TRUE", {
  tmp <- withr::local_tempdir()

  local_mocked_bindings(
    fileDownload = function(url, destfile, ...) {
      saveRDS("downloaded", destfile)
      invisible(destfile)
    }
  )

  msgs <- capture_messages(LoadRefData("survival", path = tmp, verbose = TRUE))
  txt <- paste(msgs, collapse = "\n")
  expect_match(txt, "Downloading reference data")
  expect_match(txt, "Trying GitHub Raw")
  expect_match(txt, "Data loaded successfully")
})

test_that("LoadRefData creates and uses the default cache directory", {
  root <- withr::local_tempdir()
  withr::local_envvar(HOME = root, LOCALAPPDATA = root)

  local_mocked_bindings(
    fileDownload = function(url, destfile, ...) {
      saveRDS("default-cache", destfile)
      invisible(destfile)
    }
  )

  msgs <- capture_messages(res <- LoadRefData("survival", verbose = TRUE))
  expect_match(paste(msgs, collapse = "\n"), "Using default cache path")
  expect_identical(res, "default-cache")

  cache_root <- if (.Platform$OS.type == "windows") {
    file.path(root, "SigBridgeR", "cache")
  } else {
    file.path(root, ".cache", "SigBridgeR")
  }
  expect_true(dir.exists(cache_root))
  expect_true(file.exists(file.path(cache_root, "survival_ref_data.rds")))
})
