#' @title Configure Parallel Execution Backends
#' @description
#' Unified interface to configure thread counts and acceleration features across
#' system-level (OpenMP/data.table/cheapr) and TensorFlow backends. Uses a hierarchical
#' configuration model where global \code{threads} serves as default for all
#' backends unless explicitly overridden.
#'
#' \strong{Critical}: TensorFlow configuration must be set \emph{before} importing
#' TensorFlow via \code{reticulate::import("tensorflow")}.
#'
#' @param threads Integer. Global thread count used as default for all backends.
#'   Default: `2L`
#' @param dt Integer. Thread count for data.table (default: inherited from \code{threads}).
#' @param cheapr Integer. Thread count for cheapr (default: inherited from \code{threads}).
#' @param qs2 Integer. Thread count for qs2 (default: inherited from \code{threads}).
#' @param openmp Integer. Thread count for OpenMP (default: inherited from \code{threads}).
#' @param WGCNA Integer. Thread count for WGCNA (default: inherited from \code{threads}).
#' @param tf_config Named list for TensorFlow-specific configuration:
#'   * `xla_flag`: Character. XLA JIT compilation flags (default: auto-optimized)
#'   * `xla_device`: Integer. XLA device ID (default: `1L`)
#'   * `inter_op`: Integer. Inter-op parallelism threads (default: `max(2, floor(threads/4))`)
#'   * `intra_op`: Integer. Intra-op parallelism threads (default: `1L`). If `NULL`, inherits global `threads` by default.
#' @param verbose Logical. Whether to print verbose output (default: inherited from function options).
#' @param ... Additional arguments for \code{data.table::setDTthreads()} (e.g., \code{restore}).
#'
#'
#' @return Invisible list with old/new values per backend.
#' @export
#'
#'
#' @examples
#' \dontrun{
#' # Basic usage: auto-detect and configure
#' setThreads()
#'
#' # Explicit thread count for CPU-intensive workloads
#' setThreads(threads = 12L)
#'
#' # TensorFlow-optimized configuration for deep learning
#' setThreads(
#'   threads = 8L,
#'   tf_config = list(
#'     inter_op = 2L,
#'     intra_op = 8L
#'   )
#' )
#' # Import AFTER setThreads()
#'
#' # Configure only data.table for memory-efficient workflows
#' setThreads(dt = 4L , verbose = FALSE)
#' }
setThreads <- function(
  threads = 2L,
  dt = threads,
  cheapr = threads,
  qs2 = threads,
  WGCNA = threads,
  openmp = NULL,
  tf_config = list(
    xla_flag = "--tf_xla_auto_jit=2 --tf_xla_cpu_global_jit",
    xla_device = NULL,
    # TF_NUM_INTEROP_THREADS
    inter_op = NULL, # Auto-derived: max(2, floor(threads/4))
    # TF_NUM_INTRAOP_THREADS
    intra_op = c(1L, NULL) # `NULL`: Inherits global `threads` by default
  ),
  verbose = getFuncOption("verbose"),
  ...
) {
  # ===== 1. 参数验证与默认值 =====

  chk::chk_integer(threads)
  chk::chk_list(tf_config)

  if (!is.null(tf_config$xla_flag)) {
    chk::chk_chr(tf_config$xla_flag)
  }

  if (!is.null(tf_config$inter_op)) {
    chk::chk_integer(tf_config$inter_op)
    chk::chk_range(tf_config$inter_op, c(0, Inf))
  }

  if (!is.null(tf_config$intra_op)) {
    chk::chk_integer(tf_config$intra_op)
    chk::chk_range(tf_config$intra_op, c(0, Inf))
  }

  verbose <- verbose %||% TRUE

  # ===== 2. 各后端独立配置，统一收集结果 =====
  results <- c(
    config_openmp_threads(openmp),
    config_qs2_threads(qs2),
    config_dt_threads(dt, ...),
    config_cheapr_threads(cheapr),
    config_wgcna_threads(WGCNA),
    config_tf_threads(
      threads = threads,
      tf_config = tf_config
    )
  )

  # ===== 3. 统一打印输出 =====
  if (isTRUE(verbose)) {
    print_thread_config(results)
  }

  invisible(results)
}

config_openmp_threads <- function(openmp = NULL) {
  if (!is.numeric(openmp)) {
    return(list())
  }

  old_val <- Sys.getenv("OMP_NUM_THREADS", unset = "")
  Sys.setenv(OMP_NUM_THREADS = as.character(openmp))

  list(
    openmp = list(
      name = "OMP_NUM_THREADS",
      old = if (old_val == "") "unset (1)" else old_val,
      new = openmp
    )
  )
}


config_qs2_threads <- function(qs2 = NULL) {
  if (!is.numeric(qs2)) {
    return(list())
  }

  check_installed("qs2")

  old_qs2 <- qs2::qopt("nthreads")
  qs2::qopt(parameter = "nthreads", value = qs2)

  list(
    qs2 = list(
      name = "qs2",
      old = old_qs2,
      new = qs2
    )
  )
}


config_dt_threads <- function(dt = NULL, ...) {
  if (!is.numeric(dt)) {
    return(list())
  }

  old_dt <- data.table::getDTthreads(verbose = FALSE)

  exec(
    .fn = data.table::setDTthreads,
    threads = dt,
    !!!SigBridgeRUtils::FilterArgs4Func(
      list(...),
      data.table::setDTthreads
    )
  )

  list(
    dt = list(
      name = "data.table",
      old = old_dt,
      new = dt
    )
  )
}

config_cheapr_threads <- function(cheapr = NULL) {
  if (!is.numeric(cheapr)) {
    return(list())
  }

  check_installed("cheapr")

  old_chpr <- cheapr::get_threads()
  cheapr::set_threads(cheapr)

  list(
    cheapr = list(
      name = "cheapr",
      old = if (old_chpr == 2L) "default (2)" else old_chpr,
      new = cheapr
    )
  )
}

config_wgcna_threads <- function(WGCNA = NULL) {
  if (!is.numeric(WGCNA)) {
    return(list())
  }

  check_installed("WGCNA")

  old_wgcna <- WGCNA::WGCNAnThreads()
  WGCNA::enableWGCNAThreads(WGCNA)

  list(
    WGCNA = list(
      name = "WGCNA",
      old = old_wgcna,
      new = WGCNA
    )
  )
}


config_tf_threads <- function(
  threads,
  tf_config = list(
    xla_flag = "--tf_xla_auto_jit=2 --tf_xla_cpu_global_jit",
    xla_device = NULL,
    inter_op = NULL,
    intra_op = c(1L, NULL)
  )
) {
  tf_config_default <- list(
    xla_flag = "--tf_xla_auto_jit=2 --tf_xla_cpu_global_jit",
    xla_device = NULL,
    inter_op = NULL,
    intra_op = c(1L, NULL)
  )

  tf_config <- utils::modifyList(tf_config_default, tf_config)

  tf_config$xla_device <- tf_config$xla_device %||% 1L

  tf_config$inter_op <- tf_config$inter_op %||%
    as.integer(max(2L, floor(threads / 4L)))

  tf_config$intra_op <- tf_config$intra_op %||% threads

  tf_results <- list()

  old_xla_flag <- Sys.getenv("TF_XLA_FLAGS", unset = "")
  Sys.setenv(TF_XLA_FLAGS = tf_config$xla_flag)
  tf_results$xla_flag <- list(
    name = "TF_XLA_FLAGS",
    old = if (old_xla_flag == "") "unset" else old_xla_flag,
    new = tf_config$xla_flag
  )

  old_xla_device <- Sys.getenv("TF_XLA_CPU_GLOBAL_JIT", unset = "")
  Sys.setenv(TF_XLA_CPU_GLOBAL_JIT = as.character(tf_config$xla_device))
  tf_results$xla_device <- list(
    name = "TF_XLA_CPU_GLOBAL_JIT",
    old = if (old_xla_device == "") "unset" else old_xla_device,
    new = tf_config$xla_device
  )

  old_inter <- Sys.getenv("TF_NUM_INTEROP_THREADS", unset = "")
  Sys.setenv(TF_NUM_INTEROP_THREADS = as.character(tf_config$inter_op))
  tf_results$inter_op <- list(
    name = "TF_NUM_INTEROP_THREADS",
    old = if (old_inter == "") "unset" else old_inter,
    new = tf_config$inter_op
  )

  old_intra <- Sys.getenv("TF_NUM_INTRAOP_THREADS", unset = "")
  Sys.setenv(TF_NUM_INTRAOP_THREADS = as.character(tf_config$intra_op))
  tf_results$intra_op <- list(
    name = "TF_NUM_INTRAOP_THREADS",
    old = if (old_intra == "") "unset" else old_intra,
    new = tf_config$intra_op
  )

  tf_results
}


print_thread_config <- function(results) {
  if (length(results) == 0) {
    return(invisible(NULL))
  }

  purrr::walk(results, function(cfg) {
    cli::cli_text(
      "{cli::symbol$bullet} {.field {cfg$name}}: {cfg$old} -> {cfg$new}"
    )
  })

  invisible(NULL)
}
