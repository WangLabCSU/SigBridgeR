#' @title Configure Parallel Execution Backends
#' @description
#' Unified interface to configure thread counts and acceleration features across
#' system-level (OpenMP/data.table) and TensorFlow backends. Uses a hierarchical
#' configuration model:
#' \itemize{
#'   \item \code{threads}: Global thread count (default: half physical cores)
#'   \item \code{backend}: System-level backends to configure
#'   \item \code{tf_config}: Named list for TensorFlow-specific settings
#' }
#'
#' @param threads Integer. Global thread count used as default for all backends.
#'   If \code{NULL} (default), uses \code{floor(availableCores() / 2)}.
#'   Applied to: OpenMP, data.table, and TensorFlow intra-op (unless overridden).
#' @param backend Character vector. System-level backends to configure:
#'   \code{"openmp"} (sets \code{OMP_NUM_THREADS}), \code{"dt"} (data.table threads).
#'   Default: \code{c("openmp", "dt")}.
#' @param tf_config Named list for TensorFlow-specific configuration:
#'   \itemize{
#'     \item \code{xla}: Logical. Enable XLA JIT compilation (default: \code{FALSE})
#'     \item \code{inter_op}: Integer. Inter-op parallelism threads (default: auto-derived)
#'     \item \code{intra_op}: Integer. Intra-op parallelism threads (default: \code{threads})
#'   }
#'   \strong{Must be set BEFORE importing TensorFlow}.
#' @param ... Additional arguments for \code{data.table::setDTthreads()} (e.g., \code{restore}).
#'   Includes \code{verbose} logical flag (default: \code{TRUE}).
#'
#' @return Invisible list with old/new values per backend.
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Minimal: auto-configure all backends with half cores
#' setThreads()
#'
#' # Explicit global thread count
#' setThreads(threads = 8)
#'
#' # TensorFlow optimization (MUST run BEFORE importing TensorFlow)
#' setThreads(
#'   threads = 8,
#'   tf_config = list(
#'     xla = TRUE,
#'     inter_op = 2,
#'     intra_op = 8
#'   )
#' )
#' reticulate::import("tensorflow")  # Import AFTER configuration
#'
#' # Fine-grained control: override specific backends
#' setThreads(
#'   threads = 16,
#'   backend = "dt",  # Only configure data.table
#'   tf_config = list(inter_op = 4)  # intra_op inherits threads=16
#' )
#' }
setThreads <- function(
  threads = NULL,
  backend = c("openmp", "dt"),
  tf_config = NULL,
  ...
) {
  # ===== 1. 参数验证与默认值 =====
  rlang::check_installed("future")
  backend <- tolower(backend)

  # 全局线程数（系统级默认值）
  threads <- threads %||% as.integer(floor(future::availableCores() / 2))
  chk::chk_integer(threads)

  # TensorFlow配置标准化
  tf_config <- tf_config %||% list()
  tf_config_default <- list(
    xla = FALSE,
    inter_op = NULL, # Auto-derived: max(2, floor(threads/4))
    intra_op = NULL # Inherits global `threads` by default
  )
  tf_config <- utils::modifyList(tf_config_default, tf_config)

  chk::chk_flag(tf_config$xla)
  if (!is.null(tf_config$inter_op)) {
    chk::chk_integer(tf_config$inter_op)
    chk::chk_range(tf_config$inter_op, c(0, Inf))
  }
  if (!is.null(tf_config$intra_op)) {
    chk::chk_integer(tf_config$intra_op)
    chk::chk_range(tf_config$intra_op, c(0, Inf))
  }

  # 推导未指定的TF线程数
  if (is.null(tf_config$inter_op)) {
    tf_config$inter_op <- as.integer(max(2L, floor(threads / 4)))
  }
  if (is.null(tf_config$intra_op)) {
    tf_config$intra_op <- threads
  }

  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose") %||% TRUE
  results <- rlang::list2()

  # ===== 3. 系统级后端配置 =====
  sys_results <- rlang::list2()

  if ("openmp" %in% backend) {
    old_val <- Sys.getenv("OMP_NUM_THREADS", unset = "")
    Sys.setenv(OMP_NUM_THREADS = as.character(threads))
    sys_results$openmp <- rlang::list2(
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
    sys_results$dt <- rlang::list2(
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
  tf_results <- rlang::list2()

  # XLA JIT
  if (tf_config$xla) {
    old_val <- Sys.getenv("TF_XLA_FLAGS", unset = "")
    new_val <- "--tf_xla_auto_jit=2 --tf_xla_cpu_global_jit"
    Sys.setenv(TF_XLA_FLAGS = new_val)
    tf_results$xla <- rlang::list2(
      name = "TF_XLA_FLAGS",
      old = if (old_val == "") "unset" else old_val,
      new = new_val
    )
  }

  # Inter-op threads
  old_inter <- Sys.getenv("TF_NUM_INTEROP_THREADS", unset = "")
  Sys.setenv(TF_NUM_INTEROP_THREADS = as.character(tf_config$inter_op))
  tf_results$inter_op <- rlang::list2(
    name = "TF_NUM_INTEROP_THREADS",
    old = if (old_inter == "") "unset" else old_inter,
    new = tf_config$inter_op
  )

  # Intra-op threads
  old_intra <- Sys.getenv("TF_NUM_INTRAOP_THREADS", unset = "")
  Sys.setenv(TF_NUM_INTRAOP_THREADS = as.character(tf_config$intra_op))
  tf_results$intra_op <- rlang::list2(
    name = "TF_NUM_INTRAOP_THREADS",
    old = if (old_intra == "") "unset" else old_intra,
    new = tf_config$intra_op
  )

  # 仅当用户显式请求TF配置时才显示（避免污染日志）
  has_explicit_tf_config <- !is.null(tf_config) &&
    (tf_config$xla ||
      !is.null(tf_config$inter_op) ||
      !is.null(tf_config$intra_op))

  if (has_explicit_tf_config && verbose) {
    cli::cli_h2(cli::col_yellow("TensorFlow Configuration"))
    purrr::walk(tf_results, function(cfg) {
      cli::cli_text(
        "{cli::symbol$bullet} {.field {cfg$name}}: {cfg$old} -> {cfg$new}"
      )
    })
  }

  invisible(c(results, tf_results))
}
