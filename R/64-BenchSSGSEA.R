BenchSSGSEA <- function(obj, ...) {
  UseMethod("BenchSSGSEA")
}

BenchSSGSEA.Seurat <- function(obj, ...) {}

`BenchSSGSEA.SigBridgeR::ScreenMethodResult` <- function(obj, ...) {}

BenchSSGSEA.default <- function(obj, ...) {
  cls_obj <- class(obj)
  avaialble_methods <- utils::methods(BenchSSGSEA, all.names = TRUE) |>
    as.character() |>
    gsub(pattern = ".*\\.", replacement = "", x = _)
  Abort("BenchSSGSEA is not implemented for class: {.csl {cls_obj}}", "")
}
