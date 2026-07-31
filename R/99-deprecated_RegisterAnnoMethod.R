# nocov start

#' @title Register an Annotation Method into the Strategy Registry
#'
#' @description
#' `r lifecycle::badge('deprecated')` Registers one or more user-defined annotation functions into a shared registry
#' environment (e.g., \code{SCAnnotateStrategy}). Each method is stored with minimal
#' metadata (\code{method_name} and \code{executor}) to enable dynamic dispatch in
#' annotation pipelines.
#'
#' This function supports extensible cell type annotation frameworks—ideal for integrating
#' novel classifiers (e.g., custom SingleR wrappers, marker-based assigners, or deep
#' learning models) while maintaining a uniform execution interface.
#'
#' @param ... Named arguments where each value is a function implementing an
#'   annotation method. The name becomes the method identifier (e.g.,
#'   \code{SingleR = SingleRAnnotate} registers under key \code{"SingleR"}).
#'   Unnamed arguments are auto-named using their expression via \code{get_names_4_ids()}.
#' @param registry Target registry to store annotation methods.
#'   Default: \code{SCAnnotateStrategy}.
#' @param overwrite Logical. If \code{FALSE} (default), throws an error when attempting
#'   to replace an existing method. Set to \code{TRUE} to allow updates.
#' @param verbose Logical. Whether to print a success message upon registration.
#'   Default: inherits from package option \code{getOption("SigBridgeRUtils.verbose")}.
#'
#' @return Invisibly returns \code{TRUE} on successful registration.
#' @export
#' @family Registering
#' @family Single_Cell_Annotation_Method
#'
#' @examples
#' \dontrun{
#' # Register a custom annotation wrapper
#' MyAnnotator <- function(sc, ...) {
#'   # Custom logic returning annotated Seurat/SCE object
#'   sc
#' }
#'
#' RegisterAnnoMethod(
#'   custom_annot = MyAnnotator,
#'   registry = SCAnnotateStrategy
#' )
#'
#' # Attempting to re-register without `overwrite = TRUE` will fail
#' # RegisterAnnoMethod(custom_annot = AnotherAnnotator)  # Error!
#' }
RegisterAnnoMethod <- function(
  ...,
  registry = SCAnnotateStrategy,
  overwrite = FALSE,
  verbose = getFuncOption("verbose")
) {
  lifecycle::deprecate_warn("4.0.0", "RegisterAnnoMethod()", "Register()")
  chk::chk_logical(verbose)
  chk::chk_logical(overwrite)
  chk::chk_environment(registry)

  # * detect where are functions
  dots <- rlang::list2(...)
  is_fun <- vapply(
    X = dots,
    FUN = \(x) rlang::is_function(x),
    FUN.VALUE = logical(1)
  )
  method_names <- get_names_4_ids(..., .quoses = rlang::enquos(...))[is_fun]

  for (i in seq_len(sum(is_fun))) {
    method_name <- method_names[i]
    executor <- dots[is_fun][[i]]

    # Check for existing method
    exists_already <- method_name %chin% names(registry)
    if (exists_already && !overwrite) {
      Abort(
        "Method already exists: {.val {method_name}}",
        info = "Registered methods: {.val {names(registry)}}\nUse {.code overwrite = TRUE} to force replacement"
      )
    }

    registry[[method_name]] <- rlang::list2(
      method_name = method_name,
      executor = executor
    )

    if (verbose) {
      if (exists_already) {
        cli::cli_alert_warning("Updated {.field {method_name}}")
      } else {
        cli::cli_alert_success("Registered {.field {method_name}}")
      }
    }
  }

  invisible(TRUE)
}
