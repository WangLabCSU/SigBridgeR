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
      method_version = method_version %||%
        if (is.null(py_env)) {
          r_pkg_version(pkg_name = pkg_name)
        } else {
          py_pkg_version(pkg_name = pkg_name, conda_env = py_env)
        }
    )
  }
)
