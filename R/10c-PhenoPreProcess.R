PhenoPreProcess <- new_generic(name = "PhenoPreProcess", dispatch_args = "data")

method(generic = PhenoPreProcess, class = class_datalike2d) <- function(
  data,
  ...
) {}

method(generic = PhenoPreProcess, class = class_any) <- function(data, ...) {
  cls_data <- class(Data)
  cli::cli_abort(c("x" = "Unsupported data class {..cls {cls_data}}"))
}
