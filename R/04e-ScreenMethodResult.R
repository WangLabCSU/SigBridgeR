# --------------------------------------------------------------------------------

#' @title ScreenMethodResult: Screening Method Result Container
#'
#' @description
#' `ScreenMethodResult` stores the results of a phenotype-associated cell
#' screening analysis. It wraps a Seurat object containing the screened
#' single-cell data and supports attaching additional dynamic properties
#' (e.g., screening statistics, metadata) via `...` in the constructor.
#'
#' @details
#' This class is defined using the S7 object-oriented system and inherits from
#' [SigBridgeRBase]. The core required property is `scRNA_data`, which must
#' be a [Seurat][SeuratObject::Seurat-class] object containing the screening
#' results.
#'
#' Additional properties can be attached dynamically via the constructor's
#' `...` argument. Named arguments are set as S7 properties on the returned
#' object, allowing flexible result storage without subclassing.
#'
#' ## Validator
#'
#' The class validator enforces:
#' \itemize{
#'   \item `scRNA_data` must be present on the object.
#'   \item `scRNA_data` must be a [Seurat][SeuratObject::Seurat-class] object.
#' }
#'
#' @param ... Named arguments to set as additional S7 properties on the
#'   result object (e.g., screening statistics, method-specific outputs).
#'   The `scRNA_data` property must be provided as one of the named
#'   arguments.
#'
#' @returns A `ScreenMethodResult` S7 object.
#'
#' @family S7-Classes
#' @export
ScreenMethodResult <- new_class(
  name = "ScreenMethodResult",
  parent = SigBridgeRBase,
  properties = list(
    scRNA_data = property_seurat
  ),
  validator = \(self) {
    if (!prop_exists(self, "scRNA_data")) {
      return(cli::cli_fmt(cli::cli_text(
        "{.val scRNA_data} is required"
      )))
    }

    if (!inherits(self@scRNA_data, "Seurat")) {
      return(cli::cli_fmt(cli::cli_text(
        "{.val scRNA_data} must be a {.cls Seurat} object"
      )))
    }
  },
  constructor = \(...) {
    dots <- list2(...)

    obj <- new_object(S7_object(), sigbridger_version = get_pkg_version())
    props(obj) <- dots
    obj
  }
)
