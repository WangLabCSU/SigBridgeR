# nocov start

utils::globalVariables(c(
  # AddMisc
  "base_key",
  ".",
  "suffix",
  "max_suffix",
  "..duplicate_cols"
))

#' @inherit SigBridgeRUtils::AddMisc
#' @export
AddMisc <- SigBridgeRUtils::AddMisc

#' @inherit SigBridgeRUtils::SetupPyEnv
#' @export
SetupPyEnv <- SigBridgeRUtils::SetupPyEnv

#' @inherit SigBridgeRUtils::ListPyEnv
#' @export
ListPyEnv <- SigBridgeRUtils::ListPyEnv


#' @keywords internal
ts_cli <- SigBridgeRUtils::CreateTimeStampCliEnv()

# ? Package options

#' @inherit SigBridgeRUtils::getFuncOption
#' @export
getFuncOption <- SigBridgeRUtils::getFuncOption

#' @inherit SigBridgeRUtils::setFuncOption
#' @export
setFuncOption <- SigBridgeRUtils::setFuncOption
