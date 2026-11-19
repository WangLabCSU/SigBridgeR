#' @title BenchMemoryConfig: Time and Memory Benchmark Configuration
#'
#' @description
#' `BenchMemoryConfig` stores the configuration of a time and memory benchmark
#' run performed by [BenchMemory()]. It extends [BenchConfig] with the
#' benchmarking controls collected in `bench_config`; the inputs shared by every
#' screening method (`matched_bulk`, `sc_data`, `phenotype`, `label_type`,
#' `phenotype_class`, `screen_method`) and the method specific arguments
#' collected in `method_param` are inherited from [BenchConfig].
#'
#'
#' @inheritParams BenchConfig
#' @param bench_config Named `list` of benchmarking controls (rounds,
#'   iterations, order, budgets, ...) consumed by [BenchMemory()]. Default:
#'   `list()`, meaning all defaults. Unknown fields abort.
#' @param ... Must be empty. Method specific arguments belong in `method_param`;
#'   unlike the other benchmark configurations nothing is stored in `param`.
#'
#' @returns A `BenchMemoryConfig` S7 object.
#'
#' @family S7-Classes
#' @export
BenchMemoryConfig <- S7::new_class(
  name = "BenchMemoryConfig",
  parent = BenchConfig,
  properties = list(
    bench_config = property_list
  ),
  constructor = \(
    matched_bulk,
    sc_data,
    phenotype,
    label_type = "{method_name}",
    phenotype_class = c("binary", "survival", "continuous"),
    screen_method = names(ScreenStrategy),
    method_param = MethodParamSet(),
    bench_config = BenchMemoryConfigSet(),
    ...
  ) {
    rlang::check_dots_empty()

    S7::new_object(
      S7::S7_object(),
      matched_bulk = matched_bulk,
      sc_data = sc_data,
      phenotype = phenotype,
      label_type = label_type,
      phenotype_class = phenotype_class,
      screen_method = screen_method,
      method_param = method_param,
      bench_config = bench_config,
      sigbridger_version = get_pkg_version()
    )
  }
)


#' @title Resolve Benchmark Controls
#' @description
#' Merges the user supplied `bench_config` with the defaults and validates the
#' names, so that a typo aborts instead of being silently ignored.
#' @param rounds `integer`. Number of benchmark rounds; every method is screened
#'   once per round and per-method summaries are taken across rounds. Odd values
#'   let the median ignore a single cold round. Must be at least 1. Default: `3L`.
#' @param order `character`. Order in which methods run within a round:
#'   `"random"` shuffles them (seeded by `seed`) to spread cold-cache and
#'   package-caching effects; `"fixed"` keeps the given method order. Default:
#'   `c("random", "fixed")` (first value used).
#' @param seed `integer`. Base seed for shuffling method order each round; the
#'   round number is added to it. Also set inside each measurement child process.
#'   Default: `123L`.
#' @param measure_rss `logical`. Whether to record peak process resident set size
#'   (RSS) in addition to R heap memory. Supported on Linux only. Default: `FALSE`.
#' @param on_failure `character`. Action when a method times out or errors in a
#'   round: `"skip"` warns and continues with the remaining methods, `"stop"`
#'   aborts the whole benchmark. Default: `c("skip", "stop")` (first value used).
#' @param set_threads `integer`. Number of threads the measurement child process
#'   is pinned to via [setThreads()], so runs stay comparable across methods.
#'   Default: `1L`.
#' @param timeout_sec `numeric`. Maximum wall-clock seconds allowed for one
#'   method in one round; the child process is killed when the limit is exceeded.
#'   Use `Inf` to disable the limit. Default: `10800L` (3 hours).
#' @inheritParams bench::mark
#'
#' @return A named `list` of resolved controls.
#'
#' @examples
#' \donttest{
#' BenchMemoryConfigSet()
#' }
#'
#' @family ConfigSet
#'
#' @export
BenchMemoryConfigSet <- function(
  rounds = 3L,
  min_iterations = 1L,
  max_iterations = 1L,
  order = c("random", "fixed"),
  seed = getFuncOption("seed"), # 123L
  memory = TRUE,
  filter_gc = FALSE,
  check = FALSE,
  measure_rss = FALSE,
  on_failure = c("skip", "stop"),
  set_threads = 1L,
  min_time = 0.5,
  relative = FALSE,
  time_unit = NULL,
  timeout_sec = 10800L, # 3h
  ...
) {
  arg_names <- rlang::fn_fmls_names()
  # cpp
  dots <- resolve_param_case(
    dots = list(...),
    canonical = arg_names,
    env = rlang::current_env(),
    lower_inner = FALSE
  )

  unknown <- setdiff(names(dots), arg_names)
  if (length(unknown) > 0L) {
    Abort(
      "Unknown {.arg bench_config} field(s): {.val {unknown}}",
      tips = "Valid fields: {.val {toString(arg_names)}}"
    )
  }

  rounds <- as.integer(rounds)
  min_iterations <- as.integer(min_iterations)
  max_iterations <- as.integer(max_iterations)
  order <- SigBridgeRUtils::MatchArg(order, c("random", "fixed"))
  on_failure <- SigBridgeRUtils::MatchArg(
    on_failure,
    c("skip", "stop")
  )

  if (rounds < 1L || min_iterations < 1L) {
    Abort("{.val rounds} and {.val min_iterations} must be at least 1")
  }
  seed <- as.integer(seed)

  rlang::check_bool(memory)
  rlang::check_bool(filter_gc)
  rlang::check_bool(check)
  rlang::check_bool(measure_rss)
  rlang::check_bool(relative)
  chk::chk_integer(set_threads)
  chk::chk_numeric(min_time)
  chk::chk_numeric(timeout_sec)

  list(
    rounds = rounds,
    min_iterations = min_iterations,
    max_iterations = max_iterations,
    order = order,
    seed = seed,
    memory = memory,
    filter_gc = filter_gc,
    check = check,
    measure_rss = measure_rss,
    on_failure = on_failure,
    set_threads = set_threads,
    min_time = min_time,
    relative = relative,
    time_unit = time_unit,
    timeout_sec = timeout_sec
  )
}
