# --------------------------------------------------------------------------------

#' @title BenchGeneScoreConfig: Gene Set Scoring Benchmark Configuration
#'
#' @description
#' `BenchGeneScoreConfig` stores the configuration of a gene set scoring
#' benchmark run performed by [BenchGeneScore()]. It extends [BenchConfig] with
#' the scoring controls collected in `bench_config`; the inputs shared by every
#' screening method (`matched_bulk`, `sc_data`, `phenotype`, `label_type`,
#' `phenotype_class`, `screen_method`) and the method specific arguments
#' collected in `method_param` are inherited from [BenchConfig].
#'
#' @details
#' For a `BenchGeneScoreConfig` input, [BenchGeneScore()] screens the data with
#' each method in `screen_method` and scores the accumulated cell labels in a
#' single pass, so the scores are comparable across screening methods by
#' construction. `bench_config` selects the scoring backend and collects the
#' arguments forwarded to it, as built by [BenchGeneScoreConfigSet()].
#'
#' @inheritParams BenchConfig
#' @param bench_config Named `list` of gene set scoring controls returned by
#'   [BenchGeneScoreConfigSet()]. It holds `method` (the scoring backend, one of
#'   `"ssGSEA"`, `"AUCell"` or `"UCell"`) plus any additional arguments passed
#'   to that backend. Default: `BenchGeneScoreConfigSet()`, i.e. `ssGSEA`.
#' @param ... Must be empty. Method specific arguments belong in `method_param`.
#'
#' @returns A `BenchGeneScoreConfig` S7 object.
#'
#' @examples
#' \dontrun{
#' BenchGeneScoreConfig(
#'   matched_bulk = bulk_matrix,
#'   sc_data = seurat,
#'   phenotype = phenotype,
#'   label_type = "{method_name}",
#'   phenotype_class = "binary",
#'   screen_method = "ssGSEA",
#'   method_param = list(
#'     nperm = 1000
#'   )
#' )
#' }
#'
#'
#' @family S7-Classes
#' @export
BenchGeneScoreConfig <- S7::new_class(
  name = "BenchGeneScoreConfig",
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
    bench_config = BenchGeneScoreConfigSet(),
    ...
  ) {
    rlang::check_dots_empty()

    S7::new_object(
      BenchConfig(
        matched_bulk = matched_bulk,
        sc_data = sc_data,
        phenotype = phenotype,
        label_type = label_type,
        phenotype_class = phenotype_class,
        screen_method = screen_method,
        method_param = method_param,
        sigbridger_version = get_pkg_version()
      ),
      bench_config = bench_config
    )
  }
)

#' @title Build Gene Set Scoring Benchmark Controls
#'
#' @description
#' Collects the arguments of a gene set scoring benchmark into the
#' `bench_config` list of [BenchGeneScoreConfig()]. `method` selects the scoring
#' backend, and additional named arguments are forwarded to that backend when
#' [BenchGeneScore()] computes the scores.
#'
#' @param method Character. Scoring backend to benchmark, one of `"ssGSEA"`,
#'   `"AUCell"`, `"UCell"`, `"AddModuleScore"`; The package backing the selected backend is checked with
#'   [rlang::check_installed()]. Default: `"ssGSEA"`.
#' @param ... Additional named arguments forwarded to the selected scoring
#'   backend. They are matched against the backend's formal arguments and, when
#'   a method needs two calls, split over its stages (e.g.
#'   [GSVA::ssgseaParam()] and [GSVA::gsva()] for `"ssGSEA"`). Case variants of
#'   the function's own arguments are resolved to their canonical names.
#'
#' @return A named `list` holding `method` and the extra backend arguments, to
#'   be used as the `bench_config` of [BenchGeneScoreConfig()].
#'
#' @examplesIf rlang::is_installed("GSVA")
#' BenchGeneScoreConfigSet(method = "ssGSEA")
#'
#'
#' @seealso [BenchGeneScoreConfig()], [BenchGeneScore()]
#'
#' @family ConfigSet
#' @export
BenchGeneScoreConfigSet <- function(
  method = c(
    "ssGSEA",
    "AUCell",
    "UCell",
    "AddModuleScore"
  ),
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
  backend_pkg <- switch(
    method,
    "ssGSEA" = "GSVA",
    "AUCell" = "AUCell",
    "UCell" = "UCell",
    "AddModuleScore" = "Seurat"
  )
  rlang::check_installed(backend_pkg)
  c(
    list(method = method),
    dots # pass to specific method
  )
}
