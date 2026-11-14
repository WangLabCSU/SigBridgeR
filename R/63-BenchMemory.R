#' @export
BenchMemory <- function(
  config,
  verbose = getFuncOptiions("verbose"),
  ...
) {
  rlang::check_installed(c("bench"))
  rlang::check_dots_empty()
  threads_used <- getFuncOption("nthreads")
  if (threads_used > 1L) {
    if (rlang::is_interactive()) {
      choice <- utils::askYesNo(cli::format_inline(
        "Benchmarking with multiple threads (threads = {.val {threads_used}}). Do you want to continue?"
      ))
      if (!isTRUE(choice)) {
        stop("Benchmarking aborted.")
      }
    } else {
      Abort("Benchmarking aborted due to multiple threads being used.")
    }
  }

  on.exit(gc(verbose = FALSE))
}

BenchMemory.bench <- function(
  method_param, # method_param = lots of `method_config`s
  min_time = 0.5,
  iterations = NULL,
  min_iterations = 1,
  max_iterations = 10000,
  check = TRUE,
  memory = capabilities("profmem"),
  filter_gc = TRUE,
  relative = FALSE,
  time_unit = NULL,
  exprs = NULL,
  env = parent.frame()
) {
  # * bench::mark measures wall time + memory (allocated bytes) in one pass
  bench_res <- bench::mark(
    run = {},
    min_time = 0.5,
    iterations = NULL,
    min_iterations = 1,
    max_iterations = 10000,
    check = TRUE,
    memory = capabilities("profmem"),
    filter_gc = TRUE,
    relative = FALSE,
    time_unit = NULL,
    exprs = NULL,
    env = parent.frame()
  )
}
