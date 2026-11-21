# nocov start

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
  S7::prop(object = object, name = name)
}

#' @rdname Property-visitor
#' @export
`[[.SigBridgeR::SigBridgeRBase` <- function(object, name) {
  if (missing(name)) {
    if (S7::S7_inherits(object, ScreenMethodResult)) {
      return(object@scRNA_data@meta.data)
    }
  }
  S7::prop(object = object, name = name)
}

#' Merge two objects with `+`
#'
#' `x + y` combines two screening results into a single Seurat object via
#' [MergeResult()]. Both operands may be [ScreenMethodResult] objects or Seurat
#' objects; a [ScreenMethodResult] is unwrapped to its `scRNA_data` slot before
#' merging.
#'
#'
#'
#' @param x A ScreenMethodResult or Seurat object as the main receiver.
#' @param y A ScreenMethodResult or Seurat object as the main donor.
#' @return A merged Seurat object.
#' @name plus-methods
NULL

#' @rdname plus-methods
#' @rawNamespace S3method("+","SigBridgeR::ScreenMethodResult")
#' @export
`+.SigBridgeR::ScreenMethodResult` <- function(x, y) {
  if (!S7::S7_inherits(y, ScreenMethodResult)) {
    Abort(
      "`y` must be a {.cls ScreenMethodResult} object",
      type = "[TYPE ERROR]"
    )
  }
  MergeResult(x@scRNA_data, y@scRNA_data)
}

#' @rdname plus-methods
#' @rawNamespace S3method("+",Seurat)
#' @export
`+.Seurat` <- function(x, y) {
  if (!inherits(y, "Seurat")) {
    Abort(
      "`y` must be a {.cls Seurat} object",
      type = "[TYPE ERROR]"
    )
  }
  MergeResult(x, y)
}

# nocov end
