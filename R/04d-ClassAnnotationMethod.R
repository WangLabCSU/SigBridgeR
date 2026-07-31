# --------------------------------------------------------------------------------

#' @title AnnotationMethod: Annotation Method Class for Cell Type Identification
#'
#' @description
#' `AnnotationMethod` is the S7 class for cell type annotation methods. It
#' extends [Method] with Python environment support, enabling annotation
#' methods that depend on Python-based tools (e.g., SingleR, scArches).
#'
#' @details
#' This class inherits from [Method] and adds a `py_env` property for
#' specifying a conda environment or Python path. When `py_env` is `NULL`,
#' the method version is resolved from an R package via
#' `r_pkg_version()`. When a Python environment is specified, the version
#' is resolved via `py_pkg_version()`.
#'
#' @inheritParams Method
#' @param py_env `character` or `NULL`. A conda environment name or Python
#'   path for Python-based annotation tools. When `NULL` (default), the
#'   method version is resolved from an R package.
#'
#' @returns An `AnnotationMethod` S7 object.
#'
#' @family S7-Classes
#' @export
AnnotationMethod <- new_class(
  name = "AnnotationMethod",
  properties = list(
    py_env = property_py_env
  ),
  parent = Method,
  constructor = function(
    method_name = character(),
    executor = \() ...,
    pkg_name = character(),
    method_version = character(),
    py_env = NULL
  ) {
    new_object(
      S7_object(),
      method_name = method_name,
      executor = executor,
      py_env = py_env,
      sigbridger_version = get_pkg_version(),
      method_version = method_version %||%
        if (is.null(py_env)) {
          r_pkg_version(pkg_name = pkg_name)
        } else {
          py_pkg_version(pkg_name = pkg_name, conda_env = py_env)
        }
    )
  }
)
