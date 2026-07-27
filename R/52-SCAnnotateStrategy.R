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
    CellTypist = AnnotationMethod(
      method_name = "CellTypist",
      executor = CellTypistAnnotate,
      pkg_name = "celltypist",
      method_version = "1.7.1"
    ),
    mLLMCelltype = AnnotationMethod(
      method_name = "mLLMCelltype",
      method_version = r_pkg_version("mLLMCelltype"),
      executor = mLLMCelltypeAnnotate
    ),
    SingleR = AnnotationMethod(
      method_name = "SingleR",
      method_version = r_pkg_version("SingleR"),
      executor = SingleRAnnotate
    )
  )
)
