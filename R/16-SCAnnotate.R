#' @keywords internal
SCAnnotate <- function(sc, method = c("CellTypist", "SingleR")) {
  method <- tolower(method)
  method <- SigBridgeRUtils::MatchArg(method, c("celltypist", "singler"), NULL) # must chosen
}
