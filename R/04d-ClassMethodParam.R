MethodParam <- new_class(
  name = "MethodParam",
  parent = SigBridgeRBase,
  properties = list(
    param_list = class_list,
    algorithm_version = class_character
  ),
  abstract = TRUE
)

# ------------------------------------------------------------------------------
# ! Generate Code

#' @keywords internal
GenerateMethodParamCode <- function(
  method_names = names(ScreenStrategyEnv)
) {
  func <- expression({
    METHODParam <- new_class(
      name = "METHODParam",
      properties = list(
        param_list = class_list,
        algorithm_version = class_character
      ),
      abstract = TRUE
    )
  })
  purrr::walk(
    method_names,
    ~ cli::cli_code(gsub(
      pattern = "METHOD",
      replacement = .x,
      x = as.character(func)
    ))
  )
}
