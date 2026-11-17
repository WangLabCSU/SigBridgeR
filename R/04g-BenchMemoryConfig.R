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
    screen_method = AVAILABLE_SCREEN_METHODS,
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
#' @param
#' @inheritParams bench::mark
#'
#' @return A named `list` of resolved controls.
#'
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
