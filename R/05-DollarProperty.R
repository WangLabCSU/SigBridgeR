#' Retrieve a property from SigBridgeR object
#'
#' @param object An object inheriting from SigBridgeRBase
#' @param name The name of the parameter as a character. Partial matching is not performed.
#'
#' @name Property-visitor
#' @export
NULL

#' @rdname Property-visitor
#' @export
`$.SigBridgeRBase` <- function(object, name) {
  prop(object = object, name = name)
}

#' @rdname Property-visitor
#' @export
`[[.SigBridgeRBase` <- function(object, name) {
  prop(object = object, name = name)
}
