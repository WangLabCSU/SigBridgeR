#' @title Register an Annotation Method into the Strategy Registry
#'
#' @description
#' Registers one or more user-defined annotation functions into a shared registry
#' environment (e.g., \code{SCAnnotateStrategy}). Each method is stored with minimal
#' metadata (\code{method_name} and \code{executor}) to enable dynamic dispatch in
#' annotation pipelines.
#'
#' This function supports extensible cell type annotation frameworks—ideal for integrating
#' novel classifiers (e.g., custom SingleR wrappers, marker-based assigners, or deep
#' learning models) while maintaining a uniform execution interface.
#'
#' @param ... Named arguments where each value is a **function** implementing an
#'   annotation method. The name becomes the method identifier (e.g.,
#'   \code{SingleR = SingleRAnnotate} registers under key \code{"SingleR"}).
#'   Unnamed arguments are auto-named using their expression via \code{get_names_4_ids()}.
#' @param registry Target registry to store annotation methods.
#'   Default: \code{SCAnnotateStrategy}.
#' @param verbose Logical. Whether to print a success message upon registration.
#'   Default: inherits from package option \code{getOption("SigBridgeRUtils.verbose")}.
#'
#' @return Invisibly returns \code{TRUE} on successful registration.
#' @export
#' @keywords Registering
#' @keywords Single_Cell_Annotation_Method
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
#' # Verify registration
#' names(SCAnnotateStrategy)
#' }
#'
RegisterAnnoMethod <- function(
  ...,
  registry = SCAnnotateStrategy,
  verbose = getFuncOption("verbose")
) {
  chk::chk_logical(verbose)
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

    registry[[method_name]] <- rlang::list2(
      method_name = method_name,
      executor = executor
    )
  }

  if (verbose) {
    cli::cli_alert_success(
      "Registered {.arg {method_names}}"
    )
  }

  invisible(TRUE)
}
