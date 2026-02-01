#' @keywords internal
SCAnnotate <- function(
  sc,
  method = c("CellTypist", "SingleR", "mLLMCelltype"),
  ...
) {
  chk::chk_is(sc, "Seurat")
  method <- tolower(method)
  method <- SigBridgeRUtils::MatchArg(
    method,
    c("celltypist", "singler", "mllmCelltype"),
    NULL
  ) # must chosen

  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose")
  seed <- dots$verbose %||% getFuncOption("seed")

  set.seed(seed)
}
