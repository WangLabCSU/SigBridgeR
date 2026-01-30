#' @inherit SigBridgeRUtils::AddMisc
#' @export
AddMisc <- SigBridgeRUtils::AddMisc

#' @inherit DEGAS::SetupPyEnv
#' @export
SetupPyEnv <- DEGAS::SetupPyEnv

#' @inherit DEGAS::ListPyEnv
#' @export
ListPyEnv <- DEGAS::ListPyEnv

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
