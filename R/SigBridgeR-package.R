# nocov start

## usethis namespace: start
#' @importFrom data.table := %chin% .N .SD .I
#' @importFrom Rcpp sourceCpp
#' @importFrom SigBridgeRUtils %<-%
#' @rawNamespace if (getRversion() < "4.3.0") importFrom("S7", "@")
#' @useDynLib SigBridgeR, .registration = TRUE
## usethis namespace: end
NULL

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

cli_symbol <- cli::symbol

# nocov end
