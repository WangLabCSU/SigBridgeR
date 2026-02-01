#' @title Register a Seurat Processing Strategy
#'
#' @description
#' Dynamically registers one or more new preprocessing strategies into a strategy
#' environment (i.e., \code{SCPreProcessStrategy}). Each strategy is stored under
#' a user-defined key and wraps a Seurat-compatible function to conform to the
#' standard interface: \code{function(object, params)}.
#'
#' Supports both direct function objects and character specifications of
#' package-qualified function names (e.g., \code{h = "Seurat::RunHarmony"}).
#' This enables safe, lazy resolution of external dependencies and simplifies
#' pipeline configuration via strings.
#'
#' @param ... Named arguments where each name is a **strategy key** (e.g., \code{"h"})
#'   and each value is either:
#'   \itemize{
#'     \item A **function**, or
#'     \item A **character string** specifying a function (e.g., \code{"Seurat::RunHarmony"} or \code{"stats::prcomp"}).
#'   }
#'   The function will be automatically wrapped to accept \code{object} and \code{params}, making it compatible with `SCPreProcess()` pipelines.
#' @param registry Environment. The target registry to store strategies.
#'   Default: \code{SCPreProcessStrategy}.
#' @param overwrite Logical. If \code{FALSE} (default), throws an error when attempting
#'   to replace an existing strategy. Set to \code{TRUE} to allow updates.
#' @param verbose Logical. Whether to print registration messages.
#'   Default: inherits from \code{getOption("SigBridgeR.verbose")}.
#'
#' @return Invisibly returns \code{TRUE} on successful registration.
#' @export
#'
#' @examples
#' \dontrun{
#' # Register using a package-qualified name
#' RegisterSeuratMethod(
#'   h = "Seurat::RunHarmony",
#'   registry = SCPreProcessStrategy
#' )
#'
#' # Attempting to re-register without `overwrite = TRUE` will fail
#' # RegisterSeuratMethod(h = Seurat::RunHarmony)  # Error!
#' }
#'
RegisterSeuratMethod <- function(
  ...,
  registry = SCPreProcessStrategy,
  overwrite = FALSE,
  verbose = getFuncOption("verbose")
) {
  chk::chk_logical(verbose)
  chk::chk_logical(overwrite)
  chk::chk_environment(registry)
  dots <- rlang::list2(...)
  chk::chk_named(dots)

  method_names <- names(dots)

  for (i in seq_len(length(dots))) {
    letter <- method_names[i]

    lookup <- exists(x = letter, envir = registry, inherits = FALSE)
    if (lookup && !overwrite) {
      cli::cli_abort(c(
        "x" = "Method already exists: {.val {letter}}",
        ">" = "Use `overwrite = TRUE` to force replacement"
      ))
    }
    executor <- dots[[i]]

    executor <- if (is.character(executor)) {
      parts <- strsplit(executor, "::", fixed = TRUE)[[1]]
      if (length(parts) == 2) {
        pkg <- parts[1]
        fun_name <- parts[2]
        rlang::check_installed(pkg)
        getExportedValue(pkg, fun_name)
      } else {
        match.fun(parts)
      }
    } else if (is.function(executor)) {
      executor
    } else {
      cli::cli_abort(c(
        "Function provided must be a function object or character function name."
      ))
    }

    registry[[letter]] <- function(object, params) {
      rlang::exec(executor, object = object, !!!params)
    }

    if (verbose && !lookup) {
      cli::cli_alert_success("Strategy {.field {letter}} registered")
    } else if (verbose) {
      cli::cli_alert_success("Strategy {.field {letter}} updated")
    }
  }

  invisible(TRUE)
}
