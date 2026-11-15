#' @title BenchMemoryConfig: Time and Memory Benchmark Configuration
#'
#' @description
#' `BenchMemoryConfig` stores the configuration of a time and memory benchmark
#' run performed by [BenchMemory()]. It extends [BenchConfig] with the inputs
#' shared by every screening method (`matched_bulk`, `sc_data`, `phenotype`,
#' `label_type`, `phenotype_class`, `screen_method`), the method specific
#' arguments collected in `method_param`, and the benchmarking controls
#' collected in `bench_config`.
#'
#' @details
#' Unlike the other benchmark configurations, one `BenchMemoryConfig` describes
#' several screening runs: [BenchMemory()] screens the same data with each
#' method listed in `screen_method`, one method at a time.
#'
#' * `screen_method` may only contain methods registered in [ScreenStrategy]
#'   (e.g. `"Scissor"`, `"scPAS"`, `"LP_SGL"`); the constructor aborts on any
#'   unknown name. It also aborts when the selected methods have no phenotype
#'   class in common.
#' * `method_param` is a named list with one entry per method, named
#'   `"<method>_param"` (e.g. `list(Scissor_param = list(alpha = 0.05))`).
#'   Each entry is spliced into the arguments of that method's executor, after
#'   its [ScreenMethod] mapper has run. Use [MethodParamSet()] to build the full
#'   set with empty defaults, so unspecified methods simply run with their own
#'   defaults.
#' * `label_type` is a length-one template. The default `"{method_name}"` is
#'   resolved inside [BenchMemory()] with the name of the method currently being
#'   benchmarked, so a single configuration can serve methods whose executors
#'   expect differently named phenotype labels.
#' * `bench_config` holds the controls of the measurement itself (`rounds`,
#'   `iterations`, `order`, `seed`, `warmup`, the memory/time budgets, ...).
#'   See the *Benchmark controls* section of [BenchMemory()] for the full list;
#'   unknown fields abort instead of being silently ignored.
#'
#' @param matched_bulk `matrix` or `data.frame`. Bulk expression data
#'   (genes x samples) matching `phenotype`. Shared by all methods.
#' @param sc_data A [Seurat][SeuratObject::Seurat-class] object to be screened.
#' @param phenotype Phenotype data, either a named vector or a patient survival
#'   data frame. See [BenchMemory()].
#' @param label_type `character` or `NULL`. Length-one glue-style template for
#'   the phenotype label type, with `{method_name}` available. When `NULL` or
#'   not of length one, it falls back to `"{method_name}"` with an
#'   informational message.
#' @param phenotype_class `character`. One or more of `"binary"`,
#'   `"continuous"`, `"survival"`. Default: all three.
#' @param screen_method `character`. One or more method names registered in
#'   [ScreenStrategy], giving the methods to benchmark and their order.
#'   Default: `"Scissor"`.
#' @param method_param Named `list` of method-specific arguments, each passed
#'   only to the executor of the method it is named after (e.g.
#'   `list(Scissor_param = list(alpha = 0.05))`). Default: [MethodParamSet()],
#'   i.e. one empty entry per registered method.
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
    matched_bulk = property_bench_bulk,
    sc_data = property_seurat,
    phenotype = class_any,
    label_type = property_bench_label_type,
    phenotype_class = property_phenotype_class,
    screen_method = property_chr,
    # method_param is different from method_config
    method_param = property_list,
    bench_config = property_list
  ),
  validator = \(self) {
    valid_method <- names(ScreenStrategy)

    if (!all(self@screen_method %chin% valid_method)) {
      invalid_method <- setdiff(self@screen_method, valid_method)
      return(cli::format_inline(
        "@screen_method got unknown method {.val {invalid_method}}"
      ))
    }

    # * membership is enforced by `property_phenotype_class`; several classes are
    # * allowed here because each method resolves its own supported subset
    if (length(self@phenotype_class) < 1L) {
      return(cli::format_inline(
        "@phenotype_class must contain at least one of {.val {c('binary', 'continuous', 'survival')}}"
      ))
    }

    if (length(self@method_param) < 1L) {
      return(cli::format_inline(
        "@method_param must be a {.cls non-empty list}, e.g. `list(Scissor_param = list(alpha = 0.05))`"
      ))
    }
  },
  constructor = \(
    matched_bulk,
    sc_data,
    phenotype,
    label_type = "{method_name}",
    phenotype_class = c("binary", "survival", "continuous"),
    screen_method = "Scissor",
    method_param = MethodParamSet(),
    bench_config = BenchMemoryConfigSet(),
    ...
  ) {
    rlang::check_dots_empty()
    if (is.null(label_type) || length(label_type) != 1L) {
      cli::cli_alert_info(c(
        "i" = "{.var label_type} not specified or not of length 1, using {.val {screen_method}}"
      ))
      label_type <- "{method_name}" # will be evaled in BenchMemory(), glue-style
    }

    if (!all(screen_method %chin% names(ScreenStrategy))) {
      diff <- screen_method[!screen_method %chin% names(ScreenStrategy)]
      Abort(
        "Algorithm {.val {diff}} is not supported. Please check the spelling or update the package."
      )
    }

    phenotype_classes <- SigBridgeRUtils::MatchArg(
      x = phenotype_class,
      choices = c("binary", "continuous", "survival"),
      default = NULL
    )
    phenotype_classes <- lapply(
      X = screen_method,
      FUN = function(x) ScreenStrategy[[x]]@phenotype_class
    )
    common_class <- Reduce(base::intersect, phenotype_classes)
    if (length(common_class) == 0L) {
      Abort(
        "Some of the specified algorithms do not support the phenotype class {.val {phenotype_class}}."
      )
    }

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

#' @export
MethodParamSet <- function(user_list = list()) {
  chk::chk_list(user_list)
  available_methods <- names(ScreenStrategy)
  default <- stats::setNames(
    rep(list(list()), length(available_methods)),
    paste0(available_methods, "_param")
  )
  modify_list(x = default, val = user_list) # cpp
}

#' @title Resolve Benchmark Controls
#' @description
#' Merges the user supplied `bench_config` with the defaults and validates the
#' names, so that a typo aborts instead of being silently ignored.
#' @param
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
  rlang::check_bool(auto_budget)
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
