#' @title BenchMemoryConfig: Time and Memory Benchmark Configuration
#'
#' @description
#' `BenchMemoryConfig` stores the configuration of a time and memory benchmark
#' run performed by [BenchMemory()]. It extends [BenchConfig] with the inputs
#' shared by every screening method (`matched_bulk`, `sc_data`, `phenotype`,
#' `label_type`, `phenotype_class`, `screen_method`), one parameter list per
#' screening method, and the benchmarking controls collected in `bench_config`.
#'
#' @details
#' Unlike the other benchmark configurations, one `BenchMemoryConfig` describes
#' several screening runs: [BenchMemory()] screens the same data with each
#' method listed in `screen_method`, one method at a time and in the given
#' order.
#'
#' * `screen_method` may only contain methods registered in [ScreenStrategy]
#'   (e.g. `"Scissor"`, `"scPAS"`, `"LP_SGL"`); the constructor aborts on any
#'   unknown name. It also aborts when the selected methods have no phenotype
#'   class in common.
#' * Each method owns its `<method>_param` list, which is forwarded verbatim to
#'   that method's executor (after its [ScreenMethod] mapper, if any).
#'   Methods not listed in `screen_method` simply keep the empty default.
#' * `<method>_param` names are matched case-insensitively, so `scpp_param`,
#'   `ScPP_param` and `scPP_param` are equivalent. Arguments captured by `...`
#'   that match no `<method>_param` name are reported by a warning and dropped;
#'   unlike the other benchmark configurations nothing is stored in `param`.
#' * `label_type` is a length-one glue-style template. The default
#'   `"{method_name}"` is resolved inside [BenchMemory()] with the name of the
#'   method currently being benchmarked, so a single configuration can serve
#'   methods whose executors expect differently named phenotype labels.
#' * `bench_config` holds the [bench::mark()] controls (`iterations`,
#'   `min_iterations`, `check`, ...) applied to every method.
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
#' @param method_param Named `list` of method-specific arguments, each passed only to
#'   executor of the method it is named after (e.g.
#'   `list(Scissor_param = list(alpha = 0.05))`. Default: empty list for every
#'   method.
#' @param bench_config Named `list` of benchmarking controls forwarded to
#'   [bench::mark()] by [BenchMemory()]. Default: `list()`.
#' @param ... Additional named arguments matched case-insensitively against the
#'   `<method>_param` arguments and merged into them. Unmatched names trigger a
#'   warning and are ignored.
#'
#' @returns A `BenchMemoryConfig` S7 object.
#'
#' @family S7-Classes
#' @export
BenchMemoryConfig <- new_class(
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

    if (!length(self@phenotype_class) == 1L || !all(self@phenotype_class %in% c("binary", "continuous", "survival"))) {
      return(cli::format_inline(
        "@phenotype_class must be one of {.val {c('binary', 'continuous', 'survival')}}"
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
    bench_config = list(),
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

    new_object(
      S7_object(),
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
  available_methods <- names(ScreenStrategy)
  default <- stats::setNames(
    rep(list(list()), length(available_methods)),
    paste0(available_methods, "_param")
  )
  modify_list(x = default, val = user_list)
}
