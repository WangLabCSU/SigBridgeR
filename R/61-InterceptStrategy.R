#' @title Inspect Registered Strategy Environments
#'
#' @description
#' Lists and optionally inspects the contents of predefined strategy environments
#' (e.g., \code{ScreenStrategy}, \code{SCPreProcessStrategy}, \code{AnnotationStrategy}).
#' Returns a summary of registered methods or variables, with optional detailed
#' introspection of function signatures, body structure, and closure properties.
#'
#' This utility aids in debugging, validation, and exploration of dynamically
#' registered pipelines—particularly useful when developing or auditing extensible
#' analysis frameworks.
#'
#' @param strategy Character/Environment. strategy environment to inspect.
#'   Must be one of: \code{"ScreenStrategy"}, \code{"SCPreProcessStrategy"},
#'   or \code{"SCAnnotateStrategy"}. Partial matching is supported.
#' @param ... Additional arguments (currently unused, reserved for future extension).
#'
#' @return Invisibly returns a tibble
#' @export
#'
#' @examples
#' \donttest{
#' # List all preprocessing steps
#' InterceptStrategy("SCPreProcessStrategy")
#'
#' # Get detailed function metadata for screen methods
#' details <- InterceptStrategy("ScreenStrategy")
#' }
InterceptStrategy <- function(
  strategy = c(
    "ScreenStrategy",
    "SCPreProcessStrategy",
    "SCAnnotateStrategy"
  ),
  ...
) {
  check_installed("dplyr")
  env <- if (is.character(strategy)) {
    strategy <- SigBridgeRUtils::MatchArg(
      strategy,
      c("ScreenStrategy", "SCPreProcessStrategy", "SCAnnotateStrategy"),
      NULL
    )
    get0(strategy)
  } else {
    strategy
  }

  var_names <- names(env)
  if (is.null(var_names)) {
    Abort("Found empty strategy environment: {.val {var_names}}")
  }
  info <- purrr::map(var_names, function(var_name) {
    var <- get0(strategy)[[var_name]]

    get_strategy_info(var, var_name)
  })

  dplyr::bind_rows(info)
}

get_strategy_info <- new_generic("get_strategy_info", "x")

method(get_strategy_info, SigBridgeRBase) <- function(x, name) {
  properties <- props(x)
  if (!is.null(properties$executor)) {
    properties$executor <- list(fn_fmls_names(properties$executor))
  }
  if (!is.null(properties$phenotype_class)) {
    properties$phenotype_class <- list(properties$phenotype_class)
  }
  properties$name <- name
  properties
}
method(get_strategy_info, class_function) <- function(x, name) {
  list(name = name, args = list(fn_fmls_names(x)))
}
