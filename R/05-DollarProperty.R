#' Retrieve a property from SigBridgeR object
#'
#' @param object An object inheriting from SigBridgeRBase
#' @param name The name of the parameter as a character. Partial matching is not performed.
#'
#' @name Property-visitor
NULL

#' @rdname Property-visitor
#' @export
`$.SigBridgeR::SigBridgeRBase` <- function(object, name) {
  prop(object = object, name = name)
}

#' @rdname Property-visitor
#' @export
`[[.SigBridgeR::SigBridgeRBase` <- function(object, name) {
  prop(object = object, name = name)
}

#' Merge two ScreenMethodResult objects
#' @param x A ScreenMethodResult object as main receiver
#' @param y A ScreenMethodResult object as main donor
#' @return A Seurat Object
#' @rawNamespace S3method("+","SigBridgeR::ScreenMethodResult")
#' @export
`+.SigBridgeR::ScreenMethodResult` <- function(x, y) {
  if (!S7_inherits(y, "ScreenMethodResult")) {
    Abort("`y` must be a {.cls ScreenMethodResult} object")
  }
  MergeResult(x@scRNA_data, y@scRNA_data)
}
