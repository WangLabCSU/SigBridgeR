BenchAcc <- function(obj, ...) {
  UseMethod("BenchAcc")
}

BenchAcc.Seurat <- function(obj, ...) {}

`BenchAcc.SigBridgeR::ScreenMethodResult` <- function(obj, ...) {}

BenchAcc.default <- function(obj, ...) {
  cls_obj <- class(obj)
  avaialble_methods <- utils::methods(BenchAcc, all.names = TRUE) |>
    as.character() |>
    gsub(pattern = ".*\\.", replacement = "", x = _)
  Abort("BenchAcc is not implemented for class: {.csl {cls_obj}}", "")
}
