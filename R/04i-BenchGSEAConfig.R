# --------------------------------------------------------------------------------

#' @title BenchGSEAConfig: GSEA Benchmark Configuration
#'
#' @description
#' `BenchGSEAConfig` stores the configuration of a GSEA benchmark run
#' performed by [BenchGSEA()]. It extends [BenchConfig] with the GSEA-specific
#' controls collected in `bench_config`; the inputs shared by every screening
#' method (`matched_bulk`, `sc_data`, `phenotype`, `label_type`,
#' `phenotype_class`, `screen_method`) and the method specific arguments
#' collected in `method_param` are inherited from [BenchConfig].
#'
#' @inheritParams BenchConfig
#' @param bench_config Named `list` of GSEA controls returned by
#'   [BenchGSEAConfigSet()]. It holds `group` (one group label per cell in
#'   `sc_data`) and `pathway` (a named list of gene sets), plus optional
#'   arguments forwarded to [fgsea::fgsea()] such as `minSize`, `maxSize` and
#'   `gseaParam`. `group` must have the same length as the number of cells in
#'   `sc_data`. Default: `BenchGSEAConfigSet(sc_data)`.
#' @param ... Must be empty. Method specific arguments belong in `method_param`.
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
    screen_method = names(ScreenStrategy),
    method_param = MethodParamSet(),
    bench_config = BenchGSEAConfigSet(sc_data),
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

#' @title Build GSEA Benchmark Controls
#'
#' @description
#' Collects the arguments of a GSEA benchmark into the `bench_config` list of
#' [BenchGSEAConfig()]. Arguments matching [fgsea::fgsea()] formals are handed
#' over to it; all other named arguments are kept as-is.
#'
#' @param sc_data A [Seurat][SeuratObject::Seurat-class] object to be screened.
#'   It defines the cells that `group` assigns labels to.
#' @param pathway A named `list` of gene sets, each element a character vector
#'   of gene symbols. Passed to the `pathways` argument of [fgsea::fgsea()].
#' @param ... Additional named arguments, typically further [fgsea::fgsea()]
#'   controls (e.g. `nperm`, `scoreType`). Case variants of the function's own
#'   arguments are resolved to their canonical names.
#'
#' @return A list of parameters for [fgsea::fgsea()] and [BenchGSEAConfig()].
#' @examples
#' \dontrun{
#' config <- BenchGSEAConfigSet(
#'   sc_data = sc,
#'   pathway = fgsea::examplePathways
#' )
#' }
#' @family ConfigSet
#' @export
BenchGSEAConfigSet <- function(
  sc_data, # seurat
  pathway, # pass to fgsea()
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
