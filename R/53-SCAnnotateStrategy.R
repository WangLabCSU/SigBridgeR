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
SCAnnotateStrategyEnv <- rlang::new_environment(
  list(
    CellTypist = AnnotationMethod(
      method_name = "CellTypist",
      method_version = "1.7.1",
      executor = CellTypistAnnotate,
    ),
    mLLMCelltype = AnnotationMethod(
      method_name = "mLLMCelltype",
      method_version = packageVersion("mLLMCelltype"),
      executor = mLLMCelltypeAnnotate
    ),
    SingleR = AnnotationMethod(
      method_name = "SingleR",
      method_version = packageVersion("mLLMCelltype"),
      executor = SingleRAnnotate
    )
  )
)
