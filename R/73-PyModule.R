# nocov start

#' @title Python Module Registry
#'
#' @description
#' An internal environment holding the Python modules bundled under the package's
#' `inst/python/` directory. Modules are loaded lazily at package startup via
#' [reticulate::import_from_path()] and can be accessed as `PyModule$<module>`.
#'
#' @details
#' The following Python modules are currently bundled with the package:
#'
#' `r py_module_inventory()`
#'
#' These modules provide the Python backend used by R functions such as
#' [CellTypistAnnotate()].
#'
#' @format An environment.
#' @name PyModule
#' @keywords internal
PyModule <- new.env(parent = emptyenv())

py_module_names <- function() {
  python_dir <- system.file("python", package = "SigBridgeR")

  if (!nzchar(python_dir) || !dir.exists(python_dir)) {
    return(character())
  }

  files <- list.files(
    path = python_dir,
    pattern = "[.]py$",
    recursive = TRUE,
    full.names = TRUE
  )
  files <- files[!grepl("__init__[.]py$", files)]
  if (!length(files)) {
    return(character())
  }

  modules <- substring(files, nchar(python_dir) + 2L)
  modules <- sub("[.]py$", "", modules)
  sort(gsub("/", ".", modules, fixed = TRUE))
}

py_module_inventory <- function() {
  modules <- py_module_names()
  if (!length(modules)) {
    return("* (no Python modules found)")
  }
  paste0("* `", modules, "`", collapse = "\n")
}
