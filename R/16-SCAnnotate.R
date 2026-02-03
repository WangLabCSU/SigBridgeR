#' @export
SCAnnotate <- function(
  sc,
  method = c("CellTypist", "SingleR", "mLLMCelltype"),
  ...
) {
  chk::chk_is(sc, "Seurat")
  method <- SigBridgeRUtils::MatchArg(
    method,
    names(AnnotationStrategy),
    NULL
  ) # must chosen

  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose")
  seed <- dots$verbose %||% getFuncOption("seed")

  set.seed(seed)

  AnnotationStrategy[[method]](
    sc,
    ...
  )
}
