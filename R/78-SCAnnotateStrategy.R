#' @title Registry of Cell Type Annotation Methods
#' @family Single_Cell_Annotation_Method
#' @description
#' An environment storing methods for annotating cell types.
#'
#' @details
#' Storing structure - named list
#'  - `key`: method name
#'  - `value`: list
#'     - `method_name`: method name
#'     - `executor`: function implementation of the method
#'
#' @export
SCAnnotateStrategy <- rlang::new_environment(
  list(
    CellTypist = list(
      method_name = "CellTypist",
      executor = CellTypistAnnotate
    ),
    mLLMCelltype = list(
      method_name = "mLLMCelltype",
      executor = mLLMCelltypeAnnotate
    ),
    SingleR = list(
      method_name = "SingleR",
      executor = SingleRAnnotate
    )
  )
)
