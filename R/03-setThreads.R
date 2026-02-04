#' @title Configure Parallel Thread Counts for Computational Backends
#' @description
#' Sets thread counts for OpenMP and data.table to optimize performance in
#' data-intensive workflows. Automatically uses half the available physical cores
#' if \code{n_threads} is unspecified.
#'
#' @param n_threads Integer. Number of threads to set. If \code{NULL} (default),
#'   uses \code{floor(availableCores() / 2)}.
#' @param backend Character vector. Backends to configure:
#'   \code{"openmp"} (sets \code{OMP_NUM_THREADS} env var) and/or
#'   \code{"dt"} (sets data.table threads). Default: \code{c("openmp", "dt")}.
#' @param ... Additional arguments passed to \code{data.table::setDTthreads()}
#'   (e.g., \code{restore}). Includes \code{verbose} logical flag (default: \code{TRUE}).
#'
#' @return Invisible list with old/new thread counts per backend.
#' @export
#'
#' @examples
#' \dontrun{
#' setThreads(4)                 # Set both backends to 4 threads
#' setThreads(backend = "dt")    # Auto-set only data.table threads
#' setThreads(8, verbose = FALSE) # Silent mode
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
