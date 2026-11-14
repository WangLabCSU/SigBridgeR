BenchGSEA <- function(obj, ...) {
  UseMethod("BenchGSEA")
}

BenchGSEA.Seurat <- function(obj, ...) {}

`BenchGSEA.SigBridgeR::ScreenMethodResult` <- function(obj, ...) {}

BenchGSEA.default <- function(obj, ...) {
  cls_obj <- class(obj)
  avaialble_methods <- utils::methods(BenchGSEA, all.names = TRUE) |>
    as.character() |>
    gsub(pattern = ".*\\.", replacement = "", x = _)
  Abort("BenchGSEA is not implemented for class: {.csl {cls_obj}}", "")
}
