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
#'
#' * `scRNA_data` must be present on the object.
#' * `scRNA_data` must be a [Seurat][SeuratObject::Seurat-class] object.
#'
#' @param scRNA_data A [Seurat][SeuratObject::Seurat-class] object
#'   containing the screening results. This is the core required property
#'   of the class.
#' @param ... Named arguments to set as additional S7 properties on the
#'   result object (e.g., screening statistics, method-specific outputs). Stored
#'   as a list in the `extra` property.
#' @param sigbridger_version Package Version of SigBridgeR
#'
#' @returns A `ScreenMethodResult` S7 object.
#'
#' @family S7-Classes
#' @export
ScreenMethodResult <- new_class(
  name = "ScreenMethodResult",
  parent = SigBridgeRBase,
  properties = list(
    scRNA_data = property_seurat,
    extra = property_list
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
  constructor = \(scRNA_data, ..., sigbridger_version = get_pkg_version()) {
    new_object(
      S7_object(),
      scRNA_data = scRNA_data,
      extra = list(...),
      sigbridger_version = sigbridger_version
    )
  }
)
