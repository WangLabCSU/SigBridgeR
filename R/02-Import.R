#' @inherit SigBridgeRUtils::AddMisc
#' @export
AddMisc <- SigBridgeRUtils::AddMisc

#' @inherit SigBridgeRUtils::SetupPyEnv
#' @export
SetupPyEnv <- SigBridgeRUtils::SetupPyEnv

#' @inherit SigBridgeRUtils::ListPyEnv
#' @export
ListPyEnv <- SigBridgeRUtils::ListPyEnv

# ? General global variables

#' @importFrom data.table `:=` `%chin%` `.N` `.SD`
#' @importFrom dplyr `%>%`
NULL

#' @keywords internal
ts_cli <- SigBridgeRUtils::CreateTimeStampCliEnv()

# ? Package options

#' @inherit SigBridgeRUtils::getFuncOption
#' @export
getFuncOption <- SigBridgeRUtils::getFuncOption

#' @inherit SigBridgeRUtils::setFuncOption
#' @export
setFuncOption <- SigBridgeRUtils::setFuncOption
