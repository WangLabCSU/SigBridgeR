#' @title Configure Parallel Threads
#' @description
#' Unified interface to set parallel thread counts for computational backends
#' commonly used with OpenML data workflows.
#'
#' @param n_threads Integer specifying the number of threads to use.
#'   If \code{NULL} (default), automatically use half physical cores via
#'   \code{future::availableCores()}.
#' @param backend Character vector specifying which backends to configure.
#'   Supported values:
#'   \itemize{
#'     \item{"openmp"}{ — Sets \code{OMP_NUM_THREADS}}
#'     \item{"dt"}{ — Sets \code{data.table} thread count via \code{setDTthreads()}}
#'   }
#'   Default: \code{c("openmp", "blas", "dt", "mlr3")}.
#' @param ... Additional setting, includes:
#'   \itemize{
#'     \item{\code{verbose}}{Logical. If \code{TRUE}, print messages to console.}
#'     \item{Other arguments passed to `data.table::setDTthreads()`}
#'   }
#'
#' @return Invisible list containing the actual thread counts set for each backend.
#'   Values may be integers (successfully set), \code{"env_only"} (environment
#'   variables set but no runtime control), or \code{"unsupported"} (backend not
#'   available/configurable).
#'
#' @export
#' @examples
#' \dontrun{
#' # Set all backends to 4 threads
#' setThreads(4L)
#'
#' # Configure only OpenMP backends
#' setThreads(8L, backend = c("openmp"))
#'
#' # Silent mode (no output)
#' setThreads(2L, verbose = FALSE)
#' }
setThreads <- function(
  n_threads = NULL,
  backend = c("openmp", "dt"),
  ...
) {
  rlang::check_installed("future")
  backend <- tolower(backend)
  # Determine thread count
  n_threads <- n_threads %||% floor(future::availableCores() / 2)
  chk::chk_integer(n_threads)

  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose") %||% TRUE

  if (verbose) {
    cli::cli_alert_info("Auto-detected physical cores: {.val {n_threads}}")
  }

  results <- rlang::list2()

  # 1. OpenMP environment variables
  if ("openmp" %in% backend) {
    original_OMP_NUM_THREADS <- Sys.getenv("OMP_NUM_THREADS")
    Sys.setenv(OMP_NUM_THREADS = as.character(n_threads))

    results$openmp <- rlang::list2(
      name = "OMP_NUM_THREADS",
      old = original_OMP_NUM_THREADS,
      new = n_threads
    )
  }

  # 2. data.table threads
  if ("dt" %in% backend) {
    if (verbose) {
      cli::cli_h1(cli::col_cyan("System Threads Config"))
    }
    old_dt <- data.table::getDTthreads(verbose = FALSE)
    rlang::exec(
      data.table::setDTthreads,
      threads = n_threads,
      !!!SigBridgeRUtils::FilterArgs4Func(dots, data.table::setDTthreads)
    )
    current_dt <- data.table::getDTthreads(verbose = verbose)
    results$dt <- rlang::list2(
      name = "data.table",
      old = old_dt,
      new = current_dt
    )
  }

  # Final summary
  if (verbose) {
    cli::cli_h3(cli::col_cyan("Changes"))
    purrr::walk(results, function(a_list) {
      cli::cli_text(c(
        " " = "{cli::symbol$bullet} {.val {a_list$name}}: {a_list$old} -> {a_list$new}"
      ))
    })
  }

  invisible(results)
}
