# Tests for R/63-BenchMemory_utils.R
#
# The child-process executor is always mocked: no background R session is
# spawned, so this file runs in well under a second via testthat::test_files().

.benchmem_obj <- function() {
  counts <- Matrix::Matrix(matrix(seq_len(8L), nrow = 2L), sparse = TRUE)
  suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts))
}

.benchmem_register <- function(name) {
  executor <- function(
    matched_bulk,
    sc_data,
    phenotype,
    label_type,
    phenotype_class,
    ...
  ) {
    list(scRNA_data = sc_data)
  }
  ScreenStrategy[[name]] <- ScreenMethod(
    method_name = name,
    method_version = "0.0.1",
    executor = executor,
    phenotype_class = "binary",
    mapper = NULL
  )
  withr::defer(rm(list = name, envir = ScreenStrategy))
}

.benchmem_config <- function(method, timeout_sec = 5, on_failure = "skip") {
  BenchMemoryConfig(
    matched_bulk = matrix(seq_len(4L), nrow = 2L,
                          dimnames = list(c("g1", "g2"), c("s1", "s2"))),
    sc_data = .benchmem_obj(),
    phenotype = c(c1 = 0, c2 = 1, c3 = 1, c4 = 0),
    screen_method = method,
    method_param = MethodParamSet(),
    bench_config = BenchMemoryConfigSet(
      rounds = 1L,
      order = "fixed",
      timeout_sec = timeout_sec,
      on_failure = on_failure
    )
  )
}

.benchmem_fake_bg <- function(result = NULL, die = FALSE, stay_alive = FALSE) {
  function(func, args, ...) {
    state <- new.env(parent = emptyenv())
    state$calls <- 0L
    list(
      is_alive = function() {
        state$calls <- state$calls + 1L
        stay_alive || state$calls == 1L
      },
      read_output = function() "child output",
      read_error = function() character(),
      poll_io = function(...) invisible(NULL),
      get_result = function() {
        if (die) stop("child died")
        result
      },
      kill_tree = function() invisible(NULL),
      wait = function(...) invisible(NULL)
    )
  }
}

test_that("bench_heap_used_mb returns one non-negative number", {
  used <- bench_heap_used_mb()

  expect_length(used, 1L)
  expect_type(used, "double")
  expect_gte(used, 0)
})

test_that("bench_probe_resources reports all documented fields", {
  res <- bench_probe_resources()

  expect_type(res, "list")
  expect_setequal(
    names(res),
    c("cores", "nthreads", "mem_total_gb", "mem_avail_gb", "load1", "heap_used_mb")
  )
  expect_true(is.numeric(res$load1))
  expect_true(is.numeric(res$mem_total_gb))
  expect_gte(res$heap_used_mb, 0)
})

test_that("bench_child_log keeps only the last n lines and ignores pipe errors", {
  process <- list(
    read_output = function() c("first", "second"),
    read_error = function() character()
  )

  expect_identical(bench_child_log("old", process, n = 2L), c("first", "second"))

  broken <- list(
    read_output = function() stop("gone"),
    read_error = function() stop("gone")
  )
  expect_identical(bench_child_log(c("a", "b"), broken, n = 2L), c("a", "b"))
})

test_that("fail_row builds an error row with missing memory metrics", {
  row <- fail_row("m", 2L, 3L, "timeout", "timeout", Sys.time() - 2.5)

  expect_s3_class(row, "data.table")
  expect_equal(nrow(row), 1L)
  expect_identical(row$method_name, "m")
  expect_identical(row$round, 2L)
  expect_identical(row$position, 3L)
  expect_identical(row$status, "timeout")
  expect_true(is.na(row$mem_peak_mb))
  expect_true(is.na(row$rss_peak_mb))
  expect_gte(row$elapsed, 2)
})

test_that("bench_tmp_dir stages the three shared inputs", {
  skip_if_not_installed("qs2")
  tmp_dir <- bench_tmp_dir(.benchmem_obj(), matrix(1:4, 2L), c(s1 = 0, s2 = 1))
  withr::defer(unlink(tmp_dir, recursive = TRUE, force = TRUE))

  expect_true(dir.exists(tmp_dir))
  expect_setequal(
    basename(list.files(tmp_dir)),
    c("sc_data.qs2", "matched_bulk.qdata", "phenotype.qdata")
  )
})

test_that("bench_mem_measure records a successful child measurement", {
  .benchmem_register("memok")
  config <- .benchmem_config("memok")
  result <- list(
    data.frame(
      expression = "run",
      min = 0.1,
      median = 0.1,
      stringsAsFactors = FALSE
    ),
    12.5,
    NA_real_
  )

  local_mocked_bindings(r_bg = .benchmem_fake_bg(result = result), .package = "callr")
  res <- suppressWarnings(
    bench_mem_measure(config, "memok", 1L, 1L, tempfile(), verbose = FALSE)
  )

  expect_s3_class(res, "data.table")
  expect_equal(nrow(res), 1L)
  expect_identical(res$status, "ok")
  expect_identical(res$method_name, "memok")
  expect_identical(res$expression, "memok")
  expect_identical(res$round, 1L)
  expect_identical(res$position, 1L)
  expect_equal(res$mem_peak_mb, 12.5)
  expect_true(is.na(res$error))
})

test_that("bench_mem_measure keeps a failed child as an error row", {
  .benchmem_register("memfail")
  config <- .benchmem_config("memfail")

  local_mocked_bindings(r_bg = .benchmem_fake_bg(die = TRUE), .package = "callr")
  res <- suppressWarnings(
    bench_mem_measure(config, "memfail", 1L, 1L, tempfile(), verbose = FALSE)
  )

  expect_identical(res$status, "error")
  expect_identical(res$method_name, "memfail")
  expect_match(res$error, "child process died before returning a result")
  expect_true(is.na(res$mem_peak_mb))
})

test_that("bench_mem_measure aborts on a failed child with on_failure = stop", {
  .benchmem_register("memstop")
  config <- .benchmem_config("memstop", on_failure = "stop")

  local_mocked_bindings(r_bg = .benchmem_fake_bg(die = TRUE), .package = "callr")

  expect_error(
    suppressWarnings(
      bench_mem_measure(config, "memstop", 1L, 1L, tempfile(), verbose = FALSE)
    ),
    "failed in round 1"
  )
})

test_that("bench_mem_measure records a timeout instead of hanging", {
  .benchmem_register("memtimeout")
  config <- .benchmem_config("memtimeout", timeout_sec = 0)

  local_mocked_bindings(
    r_bg = .benchmem_fake_bg(stay_alive = TRUE),
    .package = "callr"
  )
  res <- suppressWarnings(
    bench_mem_measure(config, "memtimeout", 1L, 1L, tempfile(), verbose = FALSE)
  )

  expect_identical(res$status, "timeout")
  expect_identical(res$error, "timeout")
  expect_true(is.na(res$mem_peak_mb))
})
