# setThreads ---------------------------------------------------------------

# Engines backed by optional Suggests packages are disabled by default so that
# these tests do not require qs2, cheapr or WGCNA to be installed.
call_setThreads <- function(...) {
  args <- list(
    threads = 2L,
    dt = NULL,
    qs2 = NULL,
    cheapr = NULL,
    wgcna = NULL,
    seurat = NULL,
    openmp = NULL,
    verbose = FALSE
  )
  overrides <- list(...)
  if (length(overrides) > 0L) {
    args[names(overrides)] <- overrides
  }
  do.call(setThreads, args)
}

# Every variable set by the thread helpers is reset for the duration of a test;
# withr restores the original values (including "unset") afterwards.
local_threads_env <- function(.local_envir = parent.frame()) {
  withr::local_envvar(
    OMP_NUM_THREADS = NA_character_,
    TF_XLA_FLAGS = NA_character_,
    TF_XLA_CPU_GLOBAL_JIT = NA_character_,
    TF_NUM_INTEROP_THREADS = NA_character_,
    TF_NUM_INTRAOP_THREADS = NA_character_,
    .local_envir = .local_envir
  )
}

# setThreads - input validation --------------------------------------------

test_that("setThreads rejects non-integer global thread counts", {
  expect_error(setThreads(threads = 2), class = "chk_error")
  expect_error(setThreads(threads = "2"), class = "chk_error")
  expect_error(setThreads(threads = 2.5), class = "chk_error")
})

test_that("setThreads validates the structure of tf_config", {
  expect_error(
    setThreads(threads = 2L, tf_config = "not-a-list"),
    class = "chk_error"
  )
  expect_error(
    setThreads(threads = 2L, tf_config = list(xla_flag = c("a", "b"))),
    class = "chk_error"
  )
})

test_that("setThreads validates TensorFlow thread counts", {
  expect_error(
    setThreads(threads = 2L, tf_config = list(inter_op = -1L)),
    class = "chk_error"
  )
  expect_error(
    setThreads(threads = 2L, tf_config = list(inter_op = 2)),
    class = "chk_error"
  )
  expect_error(
    setThreads(threads = 2L, tf_config = list(intra_op = -1L)),
    class = "chk_error"
  )
})

test_that("setThreads accepts valid TensorFlow overrides", {
  local_threads_env()

  expect_no_error(
    call_setThreads(tf_config = list(inter_op = 2L, intra_op = 4L))
  )
})

# setThreads - orchestration ------------------------------------------------

test_that("setThreads returns an invisible list of TensorFlow settings", {
  local_threads_env()

  out <- withVisible(call_setThreads())

  expect_false(out$visible)
  expect_named(
    out$value,
    c("xla_flag", "xla_device", "inter_op", "intra_op")
  )
  expect_identical(
    out$value$xla_flag$new,
    "--tf_xla_auto_jit=2 --tf_xla_cpu_global_jit"
  )
})

test_that("setThreads delegates each backend with its own thread count", {
  local_threads_env()
  seen <- new.env(parent = emptyenv())
  local_mocked_bindings(
    config_openmp_threads = function(openmp = NULL) {
      seen$openmp <- openmp
      list()
    },
    config_qs2_threads = function(qs2 = NULL) {
      seen$qs2 <- qs2
      list()
    },
    config_dt_threads = function(dt = NULL, ...) {
      seen$dt <- dt
      list()
    },
    config_cheapr_threads = function(cheapr = NULL) {
      seen$cheapr <- cheapr
      list()
    },
    config_wgcna_threads = function(wgcna = NULL) {
      seen$wgcna <- wgcna
      list()
    },
    config_seurat_threads = function(seurat = NULL) {
      seen$seurat <- seurat
      list()
    },
    config_tf_threads = function(threads, tf_config) {
      seen$tf <- list(threads = threads, tf_config = tf_config)
      list(tf = list(name = "TF", old = "unset", new = threads))
    },
    .package = "SigBridgeR"
  )

  res <- setThreads(threads = 6L, dt = 3L, verbose = FALSE)

  expect_equal(seen$dt, 3L)
  expect_equal(seen$qs2, 6L)
  expect_equal(seen$cheapr, 6L)
  expect_equal(seen$wgcna, 6L)
  expect_equal(seen$seurat, 6L)
  expect_null(seen$openmp)
  expect_equal(seen$tf$threads, 6L)
  expect_identical(
    seen$tf$tf_config$xla_flag,
    "--tf_xla_auto_jit=2 --tf_xla_cpu_global_jit"
  )
  expect_named(res, "tf")
})

test_that("setThreads inherits the global thread count by default", {
  local_threads_env()
  seen <- new.env(parent = emptyenv())
  local_mocked_bindings(
    config_qs2_threads = function(qs2 = NULL) {
      seen$qs2 <- qs2
      list()
    },
    config_dt_threads = function(dt = NULL, ...) {
      seen$dt <- dt
      list()
    },
    config_cheapr_threads = function(cheapr = NULL) {
      seen$cheapr <- cheapr
      list()
    },
    config_wgcna_threads = function(wgcna = NULL) {
      seen$wgcna <- wgcna
      list()
    },
    config_seurat_threads = function(seurat = NULL) {
      seen$seurat <- seurat
      list()
    },
    config_openmp_threads = function(...) list(),
    config_tf_threads = function(...) list(),
    .package = "SigBridgeR"
  )

  # Call setThreads() directly so that every backend inherits `threads`.
  setThreads(threads = 4L, verbose = FALSE)

  expect_equal(seen$dt, 4L)
  expect_equal(seen$qs2, 4L)
  expect_equal(seen$cheapr, 4L)
  expect_equal(seen$wgcna, 4L)
  expect_equal(seen$seurat, 4L)
})

test_that("setThreads forwards ... to the data.table backend", {
  local_threads_env()
  seen <- NULL
  local_mocked_bindings(
    config_dt_threads = function(dt = NULL, ...) {
      seen <<- list(dt = dt, dots = list(...))
      list()
    },
    config_qs2_threads = function(...) list(),
    config_cheapr_threads = function(...) list(),
    config_wgcna_threads = function(...) list(),
    config_seurat_threads = function(...) list(),
    config_openmp_threads = function(...) list(),
    config_tf_threads = function(...) list(),
    .package = "SigBridgeR"
  )

  call_setThreads(dt = 2L, percent = TRUE)

  expect_equal(seen$dt, 2L)
  expect_equal(seen$dots, list(percent = TRUE))
})

test_that("setThreads aggregates results from every configured backend", {
  local_threads_env()
  local_mocked_bindings(
    config_qs2_threads = function(...) list(),
    config_dt_threads = function(...) list(),
    config_cheapr_threads = function(...) list(),
    config_wgcna_threads = function(...) list(),
    config_seurat_threads = function(...) list(),
    .package = "SigBridgeR"
  )

  res <- setThreads(threads = 4L, openmp = 3L, verbose = FALSE)

  expect_named(
    res,
    c("openmp", "xla_flag", "xla_device", "inter_op", "intra_op")
  )
  expect_equal(res$openmp$new, 3L)
  expect_identical(res$inter_op$new, 2L)
})

test_that("setThreads prints only when verbose is enabled", {
  local_threads_env()
  calls <- 0L
  local_mocked_bindings(
    print_thread_config = function(results) calls <<- calls + 1L,
    .package = "SigBridgeR"
  )

  call_setThreads(verbose = FALSE)
  expect_identical(calls, 0L)

  call_setThreads(verbose = TRUE)
  expect_identical(calls, 1L)

  call_setThreads(verbose = NULL)
  expect_identical(calls, 2L)
})

test_that("setThreads passes the aggregated config to the printer", {
  local_threads_env()
  printed <- NULL
  local_mocked_bindings(
    print_thread_config = function(results) printed <<- results,
    .package = "SigBridgeR"
  )

  call_setThreads(verbose = TRUE)

  expect_named(
    printed,
    c("xla_flag", "xla_device", "inter_op", "intra_op")
  )
  expect_identical(printed$inter_op$new, 2L)
})

# config_openmp_threads ----------------------------------------------------

test_that("config_openmp_threads skips non-numeric input", {
  local_threads_env()

  expect_identical(config_openmp_threads(NULL), list())
  expect_identical(config_openmp_threads("4"), list())
  expect_identical(Sys.getenv("OMP_NUM_THREADS", unset = ""), "")
})

test_that("config_openmp_threads sets and reports the previous value", {
  local_threads_env()
  withr::local_envvar(OMP_NUM_THREADS = "7")

  res <- config_openmp_threads(4L)

  expect_identical(Sys.getenv("OMP_NUM_THREADS"), "4")
  expect_identical(res$openmp$name, "OMP_NUM_THREADS")
  expect_identical(res$openmp$old, "7")
  expect_identical(res$openmp$new, 4L)
})

test_that("config_openmp_threads reports an unset environment variable", {
  local_threads_env()

  expect_identical(config_openmp_threads(2L)$openmp$old, "unset (1)")
})

# config_dt_threads --------------------------------------------------------

test_that("config_dt_threads skips non-numeric input", {
  old <- data.table::getDTthreads(verbose = FALSE)

  expect_identical(config_dt_threads(NULL), list())
  expect_identical(config_dt_threads("2"), list())
  expect_equal(data.table::getDTthreads(verbose = FALSE), old)
})

test_that("config_dt_threads sets data.table threads and records the old value", {
  old <- data.table::getDTthreads(verbose = FALSE)
  withr::defer(data.table::setDTthreads(old))

  res <- config_dt_threads(2L)

  expect_equal(res$dt$old, old)
  expect_equal(res$dt$new, 2L)
  expect_equal(data.table::getDTthreads(verbose = FALSE), 2L)
})

test_that("config_dt_threads filters ... before calling setDTthreads", {
  seen <- NULL
  local_mocked_bindings(
    FilterArgs4Func = function(dots, .fn) {
      seen <<- dots
      list()
    },
    .package = "SigBridgeRUtils"
  )
  old <- data.table::getDTthreads(verbose = FALSE)
  withr::defer(data.table::setDTthreads(old))

  expect_no_error(config_dt_threads(2L, percent = TRUE, unsupported = 1L))

  expect_equal(seen, list(percent = TRUE, unsupported = 1L))
  expect_equal(data.table::getDTthreads(verbose = FALSE), 2L)
})

# Optional backends: qs2 / cheapr / WGCNA ----------------------------------

test_that("config_qs2_threads checks its dependency before configuring", {
  local_mocked_bindings(
    check_installed = function(pkg, ...) stop("dependency guard: ", pkg),
    .package = "SigBridgeR"
  )

  expect_identical(config_qs2_threads(NULL), list())
  expect_error(config_qs2_threads(2L), "dependency guard: qs2")
})

test_that("config_cheapr_threads checks its dependency before configuring", {
  local_mocked_bindings(
    check_installed = function(pkg, ...) stop("dependency guard: ", pkg),
    .package = "SigBridgeR"
  )

  expect_identical(config_cheapr_threads(NULL), list())
  expect_error(config_cheapr_threads(2L), "dependency guard: cheapr")
})

test_that("config_wgcna_threads checks its dependency before configuring", {
  local_mocked_bindings(
    check_installed = function(pkg, ...) stop("dependency guard: ", pkg),
    .package = "SigBridgeR"
  )

  expect_identical(config_wgcna_threads(NULL), list())
  expect_error(config_wgcna_threads(2L), "dependency guard: WGCNA")
})

# config_seurat_threads ----------------------------------------------------

test_that("config_seurat_threads skips non-numeric input", {
  expect_identical(config_seurat_threads(NULL), list())
})

test_that("config_seurat_threads warns and skips Seurat older than 6.0.0", {
  local_mocked_bindings(
    is_installed = function(...) FALSE,
    r_pkg_version = function(...) "5.2.0",
    .package = "SigBridgeR"
  )

  expect_warning(
    res <- config_seurat_threads(2L),
    "Skipping Seurat thread configuration"
  )
  expect_identical(res, list())
  expect_warning(config_seurat_threads(2L), "5.2.0")
})

test_that("config_seurat_threads configures Seurat >= 6.0.0", {
  skip_if_not(utils::packageVersion("Seurat") >= "6.0.0")

  old <- Seurat::getThreads(verbose = FALSE)
  withr::defer(Seurat::setThreads(n = old, verbose = FALSE))

  res <- config_seurat_threads(3L)

  expect_equal(res$seurat$name, "Seurat")
  expect_equal(res$seurat$old, old)
  expect_equal(res$seurat$new, 3L)
  expect_equal(Seurat::getThreads(verbose = FALSE), 3L)
})

# config_tf_threads --------------------------------------------------------

test_that("config_tf_threads applies defaults and sets environment variables", {
  local_threads_env()

  res <- config_tf_threads(8L)

  expect_named(res, c("xla_flag", "xla_device", "inter_op", "intra_op"))
  expect_identical(
    res$xla_flag$new,
    "--tf_xla_auto_jit=2 --tf_xla_cpu_global_jit"
  )
  expect_identical(res$xla_device$new, 1L)
  expect_identical(res$inter_op$new, 2L)
  expect_identical(res$intra_op$new, 1L)
  expect_identical(res$inter_op$old, "unset")
  expect_identical(
    Sys.getenv("TF_XLA_FLAGS"),
    "--tf_xla_auto_jit=2 --tf_xla_cpu_global_jit"
  )
  expect_identical(Sys.getenv("TF_XLA_CPU_GLOBAL_JIT"), "1")
  expect_identical(Sys.getenv("TF_NUM_INTEROP_THREADS"), "2")
  expect_identical(Sys.getenv("TF_NUM_INTRAOP_THREADS"), "1")
})

test_that("config_tf_threads honours explicit overrides", {
  local_threads_env()

  res <- config_tf_threads(
    8L,
    list(
      xla_flag = "--custom-xla",
      xla_device = 3L,
      inter_op = 5L,
      intra_op = 6L
    )
  )

  expect_identical(res$xla_flag$new, "--custom-xla")
  expect_identical(res$xla_device$new, 3L)
  expect_identical(res$inter_op$new, 5L)
  expect_identical(res$intra_op$new, 6L)
  expect_identical(Sys.getenv("TF_XLA_FLAGS"), "--custom-xla")
  expect_identical(Sys.getenv("TF_XLA_CPU_GLOBAL_JIT"), "3")
  expect_identical(Sys.getenv("TF_NUM_INTEROP_THREADS"), "5")
  expect_identical(Sys.getenv("TF_NUM_INTRAOP_THREADS"), "6")
})

test_that("config_tf_threads derives inter_op from the global thread count", {
  local_threads_env()

  expect_identical(config_tf_threads(2L)$inter_op$new, 2L)
  expect_identical(config_tf_threads(16L)$inter_op$new, 4L)
})

test_that("config_tf_threads inherits global threads when intra_op is NULL", {
  local_threads_env()

  res <- config_tf_threads(8L, list(intra_op = NULL))

  expect_identical(res$intra_op$new, 8L)
  expect_identical(Sys.getenv("TF_NUM_INTRAOP_THREADS"), "8")
})

test_that("config_tf_threads reports existing environment values as old", {
  local_threads_env()
  withr::local_envvar(
    TF_XLA_FLAGS = "--old-flag",
    TF_XLA_CPU_GLOBAL_JIT = "9",
    TF_NUM_INTEROP_THREADS = "3",
    TF_NUM_INTRAOP_THREADS = "4"
  )

  res <- config_tf_threads(8L)

  expect_identical(res$xla_flag$old, "--old-flag")
  expect_identical(res$xla_device$old, "9")
  expect_identical(res$inter_op$old, "3")
  expect_identical(res$intra_op$old, "4")
})

# print_thread_config ------------------------------------------------------

test_that("print_thread_config is silent and invisible for empty results", {
  calls <- 0L
  local_mocked_bindings(
    cli_li = function(...) calls <<- calls + 1L,
    .package = "cli"
  )

  out <- withVisible(print_thread_config(list()))

  expect_false(out$visible)
  expect_null(out$value)
  expect_identical(calls, 0L)
})

test_that("print_thread_config reports each backend with old and new values", {
  seen <- list()
  templates <- NULL
  local_mocked_bindings(
    cli_li = function(...) {
      templates <<- c(templates, unlist(list(...)))
      cfg <- get("cfg", envir = parent.frame(), inherits = TRUE)
      seen <<- c(seen, list(cfg))
    },
    .package = "cli"
  )

  print_thread_config(list(
    dt = list(name = "data.table", old = 8L, new = 2L),
    openmp = list(name = "OMP_NUM_THREADS", old = "unset (1)", new = 4L)
  ))

  expect_length(seen, 2L)
  expect_equal(seen[[1L]]$name, "data.table")
  expect_equal(seen[[1L]]$old, 8L)
  expect_equal(seen[[1L]]$new, 2L)
  expect_equal(seen[[2L]]$name, "OMP_NUM_THREADS")
  expect_match(templates[[1L]], "old", fixed = TRUE)
  expect_match(templates[[1L]], "new", fixed = TRUE)
})

test_that("print_thread_config returns invisibly", {
  local_mocked_bindings(cli_li = function(...) NULL, .package = "cli")

  out <- withVisible(print_thread_config(list(
    a = list(name = "A", old = "1", new = 2L)
  )))

  expect_false(out$visible)
  expect_null(out$value)
})
