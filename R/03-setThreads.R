#' @title Configure Parallel Execution Backends
#' @description
#' Unified interface to configure thread counts and acceleration features across
#' system-level (OpenMP/data.table) and TensorFlow backends. Uses a hierarchical
#' configuration model where global \code{threads} serves as default for all
#' backends unless explicitly overridden.
#'
#' \strong{Critical}: TensorFlow configuration must be set \emph{before} importing
#' TensorFlow via \code{reticulate::import("tensorflow")}.
#'
#' @param threads Integer. Global thread count used as default for all backends.
#'   If \code{NULL} (default), uses \code{floor(availableCores() / 2)}. Applied to:
#'   OpenMP, data.table, and TensorFlow intra-op (unless overridden).
#' @param backend Character vector. System-level backends to configure:
#'   \code{"openmp"} (sets \code{OMP_NUM_THREADS}), \code{"dt"} (data.table threads).
#'   Default: \code{c("openmp", "dt")}.
#' @param tf_config Named list for TensorFlow-specific configuration:
#'   \itemize{
#'     \item \code{xla_flag}: Character. XLA JIT compilation flags (default: auto-optimized)
#'     \item \code{xla_device}: Integer. XLA device ID (default: \code{1L})
#'     \item \code{inter_op}: Integer. Inter-op parallelism threads (default: \code{max(2, floor(threads/4))})
#'     \item \code{intra_op}: Integer. Intra-op parallelism threads (default: \code{1L}). If \code{NULL}, inherits global \code{threads} by default.
#'   }
#' @param ... Additional arguments for \code{data.table::setDTthreads()} (e.g., \code{restore}).
#'   Includes \code{verbose} logical flag (default: \code{TRUE}).
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
#' library(tensorflow)  # Import AFTER setThreads()
#'
#' # Configure only data.table for memory-efficient workflows
#' setThreads(threads = 4, backend = "dt", verbose = FALSE)
#' }
setThreads <- function(
  threads = NULL,
  backend = c("openmp", "dt"),
  tf_config = list(
    xla_flag = "--tf_xla_auto_jit=2 --tf_xla_cpu_global_jit",
    xla_device = NULL,
    # TF_NUM_INTEROP_THREADS
    inter_op = NULL, # Auto-derived: max(2, floor(threads/4))
    # TF_NUM_INTRAOP_THREADS
    intra_op = c(1L, NULL) # `NULL`: Inherits global `threads` by default
  ),
  ...
) {
  # ===== 1. 参数验证与默认值 =====
  rlang::check_installed("future")
  backend <- tolower(backend)

  # 全局线程数（系统级默认值）
  threads <- threads %||% as.integer(floor(future::availableCores() / 2))
  chk::chk_integer(threads)
  chk::chk_list(tf_config)

  tf_config_default <- list(
    xla_flag = "--tf_xla_auto_jit=2 --tf_xla_cpu_global_jit",
    xla_device = NULL,
    inter_op = NULL,
    intra_op = c(1L, NULL)
  )

  # TensorFlow配置标准化
  tf_config <- utils::modifyList(tf_config_default, tf_config)

  chk::chk_chr(tf_config$xla_flag)
  if (!is.null(tf_config$inter_op)) {
    chk::chk_integer(tf_config$inter_op)
    chk::chk_range(tf_config$inter_op, c(0, Inf))
  }
  if (!is.null(tf_config$intra_op)) {
    chk::chk_integer(tf_config$intra_op)
    chk::chk_range(tf_config$intra_op, c(0, Inf))
  }

  # 推导未指定的TF线程数
  tf_config$xla_device <- tf_config$xla_device %||% 1L

  tf_config$inter_op <- tf_config$inter_op %||%
    as.integer(max(2L, floor(threads / 4L)))

  tf_config$intra_op <- tf_config$intra_op %||% threads

  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose") %||% TRUE
  results <- list()

  # ===== 3. 系统级后端配置 =====
  sys_results <- list()

  if ("openmp" %in% backend) {
    old_val <- Sys.getenv("OMP_NUM_THREADS", unset = "")
    Sys.setenv(OMP_NUM_THREADS = as.character(threads))
    sys_results$openmp <- list(
      name = "OMP_NUM_THREADS",
      old = if (old_val == "") "unset" else old_val,
      new = threads
    )
  }

  if ("dt" %in% backend) {
    old_dt <- data.table::getDTthreads(verbose = FALSE)
    rlang::exec(
      data.table::setDTthreads,
      threads = threads,
      !!!SigBridgeRUtils::FilterArgs4Func(dots, data.table::setDTthreads)
    )
    new_dt <- data.table::getDTthreads(verbose = verbose)
    sys_results$dt <- list(
      name = "data.table",
      old = old_dt,
      new = new_dt
    )
  }

  if (length(sys_results) > 0 && verbose) {
    cli::cli_h2(cli::col_cyan("System-Level Configuration"))
    purrr::walk(sys_results, function(cfg) {
      cli::cli_text(
        "{cli::symbol$bullet} {.field {cfg$name}}: {cfg$old} -> {cfg$new}"
      )
    })
  }

  results <- c(results, sys_results)

  # ===== 4. TensorFlow专属配置（统一命名空间） =====
  tf_results <- list()

  # XLA JIT
  old_xla_flag <- Sys.getenv("TF_XLA_FLAGS", unset = "")
  Sys.setenv(TF_XLA_FLAGS = tf_config$xla_flag)
  tf_results$xla_flag <- list(
    name = "TF_XLA_FLAGS",
    old = if (old_val == "") "unset" else old_val,
    new = tf_config$xla_flag
  )

  old_xla_device <- Sys.getenv("TF_XLA_CPU_GLOBAL_JIT", unset = "")
  Sys.setenv(TF_XLA_CPU_GLOBAL_JIT = as.character(tf_config$xla_device))
  tf_results$xla_device <- list(
    name = "TF_XLA_CPU_GLOBAL_JIT",
    old = if (old_val == "") "unset" else old_val,
    new = tf_config$xla_device
  )

  # Inter-op threads
  old_inter <- Sys.getenv("TF_NUM_INTEROP_THREADS", unset = "")
  Sys.setenv(TF_NUM_INTEROP_THREADS = as.character(tf_config$inter_op))
  tf_results$inter_op <- list(
    name = "TF_NUM_INTEROP_THREADS",
    old = if (old_inter == "") "unset" else old_inter,
    new = tf_config$inter_op
  )

  # Intra-op threads
  old_intra <- Sys.getenv("TF_NUM_INTRAOP_THREADS", unset = "")
  Sys.setenv(TF_NUM_INTRAOP_THREADS = as.character(tf_config$intra_op))
  tf_results$intra_op <- list(
    name = "TF_NUM_INTRAOP_THREADS",
    old = if (old_intra == "") "unset" else old_intra,
    new = tf_config$intra_op
  )

  if (verbose) {
    cli::cli_h2(cli::col_yellow("TensorFlow Configuration"))
    purrr::walk(tf_results, function(cfg) {
      cli::cli_text(
        "{cli::symbol$bullet} {.field {cfg$name}}: {cfg$old} -> {cfg$new}"
      )
    })
  }

  invisible(c(results, tf_results))
}
