# --------------------------------------------------------------------------------

#' @title BenchConfig: Base Configuration for Benchmarks
#'
#' @description
#' `BenchConfig` is the common parent of the benchmark configurations in
#' SigBridgeR. It stores the inputs shared by every screening method
#' (`matched_bulk`, `sc_data`, `phenotype`, `label_type`, `phenotype_class`,
#' `screen_method`), the method specific arguments collected in `method_param`,
#' and an optional `param` list used by subclasses that forward extra arguments.
#'
#' @details
#' A single configuration may describe several screening runs: the benchmark
#' screens the same data with each method listed in `screen_method`, one method
#' at a time.
#'
#' * `screen_method` may only contain methods registered in [ScreenStrategy]
#'   (e.g. `"Scissor"`, `"scPAS"`, `"LP_SGL"`); the constructor of the concrete
#'   configuration aborts on any unknown name. It also aborts when the selected
#'   methods have no phenotype class in common.
#' * `method_param` is a named list with one entry per method, named
#'   `"<method>_param"` (e.g. `list(Scissor_param = list(alpha = 0.05))`).
#'   Each entry is spliced into the arguments of that method's executor, after
#'   its [ScreenMethod] mapper has run. Use [MethodParamSet()] to build the full
#'   set with empty defaults, so unspecified methods simply run with their own
#'   defaults.
#' * `label_type` is a length-one template. The default `"{method_name}"` is
#'   resolved by the benchmark with the name of the method currently being
#'   benchmarked, so a single configuration can serve methods whose executors
#'   expect differently named phenotype labels.
#'
#' @param matched_bulk `matrix` or `data.frame`. Bulk expression data
#'   (genes x samples) matching `phenotype`. Shared by all methods.
#' @param sc_data A [Seurat][SeuratObject::Seurat-class] object to be screened.
#' @param phenotype Phenotype data, either a named vector or a patient survival
#'   data frame. See [Screen()].
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
#' @param sigbridger_version Package Version of SigBridgeR
#' @param ... unused
#'
#' @returns A `BenchConfig` S7 object.
#'
#' @family S7-Classes
#' @export
BenchConfig <- S7::new_class(
  name = "BenchConfig",
  parent = SigBridgeRBase,
  properties = list(
    matched_bulk = property_matrix,
    sc_data = property_seurat,
    phenotype = property_phenotype,
    label_type = property_chr,
    phenotype_class = property_chr,
    screen_method = property_chr,
    # method_param is a list of method_configs
    method_param = property_data_list
  ),
  validator = \(self) {
    valid_method <- AVAILABLE_SCREEN_METHODS

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
    screen_method = AVAILABLE_SCREEN_METHODS,
    method_param = MethodParamSet(),
    sigbridger_version = get_pkg_version(),
    ...
  ) {
    rlang::check_dots_empty()
    # * check label_type
    if (is.null(label_type) || length(label_type) != 1L) {
      cli::cli_alert_info(c(
        "i" = "{.var label_type} not specified or not of length 1, using {.val {screen_method}}"
      ))
      label_type <- "{method_name}" # will be evaled in BenchMemory(), glue-style
    }

    # * check screen_method
    if (!all(screen_method %chin% AVAILABLE_SCREEN_METHODS)) {
      diff <- screen_method[!screen_method %chin% AVAILABLE_SCREEN_METHODS]
      Abort(
        "Algorithm {.val {diff}} is not supported. Please check the spelling or update the package."
      )
    }

    # * check phenotype_class
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
      sigbridger_version = sigbridger_version
    )
  }
)

#' A List of Method Parameters
#'
#' @description
#' Collects the method-specific arguments of a benchmark configuration.
#' Each method registered in [ScreenStrategy] has a corresponding argument
#' named `"<method>_param"`, which defaults to `list()`.
#'
#' @eval GenerateRoxy.MethodParamSet()
#'
#' @returns A modified list of method-specific parameters.
#'
#' @example
#' \donttest{
#' MethodParamSet(
#'   Scissor_param = list(alpha = 0.05),
#'   scPAS_param = list(alpha = 0.05)
#' )
#' }
#' @family ConfigSet
#' @export
MethodParamSet <- function(
  Scissor_param = list(),
  scAB_param = list(),
  scPAS_param = list(),
  scPP_param = list(),
  DEGAS_param = list(),
  DEGASv2_param = list(),
  LP_SGL_param = list(),
  PIPET_param = list(),
  SIDISH_param = list(),
  TiRank_param = list(),
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
      "Unknown {.arg method_param} field(s): {.val {unknown}}",
      tips = "Valid fields: {.val {toString(arg_names)}}"
    )
  }

  list(
    Scissor_param = Scissor_param,
    scAB_param = scAB_param,
    scPAS_param = scPAS_param,
    scPP_param = scPP_param,
    DEGAS_param = DEGAS_param,
    DEGASv2_param = DEGASv2_param,
    LP_SGL_param = LP_SGL_param,
    PIPET_param = PIPET_param,
    SIDISH_param = SIDISH_param,
    TiRank_param = TiRank_param
  )
}

#' Generate roxygen lines for the `@param` of [MethodParamSet()]
#'
#' One `@param` tag per registered method, in registry order.
#'
#' @returns `character()` of roxygen lines, one per tag.
#' @name GenerateRoxy
#' @keywords internal
GenerateRoxy.MethodParamSet <- function() {
  methods <- names(ScreenStrategy)
  paste0(
    "#' @param ",
    methods,
    "_param Extra arguments for method ",
    methods
  )
}
