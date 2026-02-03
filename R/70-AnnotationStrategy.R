#' @title Registry of Cell Type Annotation Methods
#' @keywords Single_Cell_Annotation_Method
#' @description
#' An object storing methods for annotating cell types.
#'
#' @details
#' Storing structure - named list
#'  - `key`: method name
#'  - `value`: list
#'     - `method_name`: method name
#'     - `executor`: function implementation of the method
#'
#' @export
AnnotationStrategy <- structure(
  list(
    CellTypist = list(
      method_name = "CellTypist",
      executor = CellTypistAnnotate
    ),
    mLLMCelltype = list(
      method_name = "mLLMCelltype",
      executor = mLLMCellTypeAnnotate
    ),
    SingleR = list(
      method_name = "SingleR",
      executor = SingleRAnnotate
    )
  ),
  class = c("AnnotationStrategy", "list")
)
