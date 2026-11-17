# --------------------------------------------------------------------------------

#' @title BenchGeneScoreConfig: ssGSEA Benchmark Configuration
#'
#' @description
#' `BenchGeneScoreConfig` stores the configuration of an ssGSEA benchmark run
#' performed by [BenchGeneScore()]. It extends [BenchConfig] with the benchmarked
#' object, and forwards every additional argument to [BenchGeneScore()] as
#' `param`.
#'
#' @inheritParams BenchGSEAConfig
#'
#' @returns A `BenchGeneScoreConfig` S7 object.
#'
#' @family S7-Classes
#' @export
BenchGeneScoreConfig <- S7::new_class(
  name = "BenchGSEAConfig",
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
    bench_config = BenchGeneScoreConfigSet(),
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
#' @param method The name of the method to benchmark.
#' @export
BenchGeneScoreConfigSet <- function(
  method = AVAILABLE_GENE_SCORE_METHODS,
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

  method <- rlang::arg_match(method)
  rlang::check_installed(method)
  c(
    list(method = method),
    dots # pass to specific method
  )
}

GENE_SCORE_METHODS <- c(
  "ssGSEA",
  "AUCell",
  "UCell"
)
