SCSCPreProcess <- new_generic(
  name = "SCPreProcess",
  dispatch_args = "data",
  fun = \(data, ...) {
    S7_dispatch()
  }
)
# 返回一个SCInput对象

method(generic = SCSCPreProcess, class = class_seurat) <- function(data, ...) {}

method(generic = SCPreProcess, class = class_InMemoryAnnData) <- function(
  data,
  ...
) {}

method(generic = SCPreProcess, class = class_AnnDataR6) <- function(
  data,
  ...
) {}

method(generic = SCPreProcess, class = class_datalike2d) <- function(
  data,
  ...
) {}

method(generic = SCPreProcess, class = class_any) <- function(data, ...) {
  cls_data <- class(Data)
  cli::cli_abort(c("x" = "Unsupported data class {..cls {cls_data}}"))
}
