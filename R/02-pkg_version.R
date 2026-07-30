#' Get Installed Version of an R Package
#'
#' @description
#' Queries the installed version of an R package via
#' [utils::packageVersion()]. Returns `NULL` on error (e.g., package not
#' installed).
#'
#' @param pkg_name A character string. The name of the R package.
#'
#' @return A character string with the version, or `NULL` if unavailable.
#'
#' @keywords internal
#' @family get-pkg-version
r_pkg_version <- function(pkg_name) {
  tryCatch(as.character(utils::packageVersion(pkg_name)), error = NULL)
}

#' Get SigBridgeR Package Version
#'
#' @description
#' Returns the installed version of the SigBridgeR package. When the package
#' is loaded via [devtools::load_all()] (i.e., in development mode), falls
#' back to reading the `Version` field from the `DESCRIPTION` file in the
#' working directory.
#'
#' @return A character string with the SigBridgeR version.
#'
#' @keywords internal
#' @family get-pkg-version
get_pkg_version <- function() {
  tryCatch(
    as.character(utils::packageVersion("SigBridgeR")),
    error = function(e) {
      # devtools::load_all() switches working directory to the package root
      read.dcf("DESCRIPTION")[, "Version"]
    }
  )
}

#' Get Installed Version of a Python Package in a Conda Environment
#'
#' @description
#' Queries the installed version of a Python package in a specified conda
#' environment by running `conda list -n <env> <pkg>`. Returns `NULL` if the
#' package is not found or conda is unavailable.
#'
#' @param pkg_name A character string. The name of the Python package.
#' @param conda_env A character string. The name of the conda environment.
#'
#' @return A character string with the version, or `NULL` if unavailable.
#'
#' @keywords internal
#' @family get-pkg-version
py_pkg_version <- function(pkg_name, conda_env) {
  # Prefer conda metadata: conda list -n ENV PKG
  conda_out <- tryCatch(
    system2(
      command = "conda",
      args = c(
        "list",
        "-n",
        shQuote(conda_env),
        shQuote(pkg_name)
      ),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character(0L)
  )

  conda_out <- conda_out[!grepl("^#", trimws(conda_out))]

  if (length(conda_out) == 0L) {
    return(NULL)
  }

  fields <- strsplit(trimws(conda_out), "\\s+") # list

  hit <- fields[vapply(
    X = fields,
    FUN = function(x) {
      length(x) >= 2L && identical(x[[1L]], pkg_name)
    },
    FUN.VALUE = logical(1L)
  )]

  if (length(hit) > 0L) {
    return(hit[[1L]][[2L]])
  }
  NULL
}
