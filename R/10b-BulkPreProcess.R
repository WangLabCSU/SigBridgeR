BulkPreProcess <- new_generic(name = "BulkPreProcess", dispatch_args = "data")

method(BulkPreProcess, class = class_datalike2d) <- function(data, ...) {}

method(BulkPreProcess, class_any) <- function(data, ...) {
  cls_data <- class(Data)
  cli::cli_abort(c("x" = "Unsupported data class {..cls {cls_data}}"))
}
