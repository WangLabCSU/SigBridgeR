# 同时用作screen result和cache config
ScreenMethodConfig <- new_class(
  name = "ScreenMethodConfig",
  parent = SigBridgeRBase,
  properties = list(
    method_name = class_character,
    method_version = class_character,
    config = class_list
  ),
  abstract = TRUE,
  validator = \(self) {
    if (!is_named(self@config)) {
      return(
        "@config must have names (sample names)"
      )
    }

    if (length(self@config) == 0L) {
      return("@config is empty")
    }
  }
)

# ------------------------------------------------------------------------------
# ! Generate Code

#' @keywords internal
GenerateMethodConfigCode <- function(
  method_names = names(ScreenStrategyEnv)
) {
  func <- expression({
    METHODConfig <- new_class(
      name = "METHODConfig",
      parent = ScreenMethodConfig
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
