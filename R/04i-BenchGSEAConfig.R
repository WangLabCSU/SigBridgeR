# --------------------------------------------------------------------------------

#' @title BenchGSEAConfig: GSEA Benchmark Configuration
#'
#' @description
#' `BenchGSEAConfig` stores the configuration of a GSEA benchmark run
#' performed by [BenchGSEA()]. It extends [BenchConfig] with the benchmarked
#' object, and forwards every additional argument to [BenchGSEA()] as `param`.
#'
#' @inheritParams BenchAccConfig
#'
#' @returns A `BenchGSEAConfig` S7 object.
#'
#' @family S7-Classes
#' @export
BenchGSEAConfig <- S7::new_class(
  name = "BenchGSEAConfig",
  parent = BenchConfig,
  properties = list(
    bench_config = property_list
  ),
  validator = \(self) {
    if (length(self@bench_config$group) != ncol(self@sc_data)) {
      Abort(
        "{.arg group} must have the same length as the number of cells",
        "Each label represents a cell"
      )
    }
  },
  constructor = \(
    matched_bulk,
    sc_data,
    phenotype,
    label_type = "{method_name}",
    phenotype_class = c("binary", "survival", "continuous"),
    screen_method = AVAILABLE_SCREEN_METHODS,
    method_param = MethodParamSet(),
    bench_config = BenchGSEAConfigSet(sc_data),
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

#'
#' @param group A vector of length number of cells, with each element being the group
#' @inheritParams fgsea::fgsea
#'
#' @return A list of parameters for `fgsea()` and `BenchGSEAConfig`
#' @example
#' \donttest{
#' BenchGSEAConfigSet()
#' }
#'
#' @export
BenchGSEAConfigSet <- function(
  sc_data, # seurat
  pathway, # pass to fgsea()
  minSize = 1,
  maxSize = length(group) - 1,
  gseaParam = 1,
  ...
) {
  rlang::check_installed("dplyr")
  arg_names <- rlang::fn_fmls_names()
  # cpp
  dots <- resolve_param_case(
    dots = list(...),
    canonical = arg_names,
    env = rlang::current_env(),
    lower_inner = FALSE
  )

  c(
    list(
      pathway = pathway # named list, elements are chr vecs
    ),
    dots # pass to fgsea()
  )
}
