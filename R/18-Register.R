#' @title Unified Registration Interface for Strategy Methods
#'
#' @description
#' A convenience wrapper that dispatches registration requests to the appropriate
#' strategy-specific registrar based on the target \code{registry}. This function
#' provides a single entry point for registering screening, preprocessing, or
#' annotation methods without needing to call \code{RegisterScreenMethod()},
#' \code{RegisterSeuratMethod()}, or \code{RegisterAnnoMethod()} directly.
#'
#' Internally routes to:
#' \itemize{
#'   \item \code{registry = "ScreenStrategy"} → \code{\link{RegisterScreenMethod}}
#'   \item \code{registry = "SCPreProcessStrategy"} → \code{\link{RegisterSeuratMethod}}
#'   \item \code{registry = "AnnotationStrategy"} → \code{\link{RegisterAnnoMethod}}
#' }
#'
#' @param ... Arguments passed to the underlying registrar. The exact requirements
#'   depend on the target \code{registry}:
#'   \describe{
#'     \item{\code{ScreenStrategy}}{Named functions with optional \code{supported_phenotypes},
#'       \code{parameter_mapper}, etc. (see \code{\link{RegisterScreenMethod}}).}
#'     \item{\code{SCPreProcessStrategy}}{Named functions or character specifications
#'       (e.g., \code{"h" = "Seurat::RunHarmony"}; see \code{\link{RegisterSeuratMethod}}).}
#'     \item{\code{AnnotationStrategy}}{Named annotation functions (see \code{\link{RegisterAnnoMethod}}).}
#'   }
#' @param registry Character. Target strategy environment for registration.
#'   Must be one of: \code{"ScreenStrategy"}, \code{"SCPreProcessStrategy"}, or
#'   \code{"AnnotationStrategy"}. Partial matching is supported (e.g., \code{"screen"} → \code{"ScreenStrategy"}).
#' @param verbose Logical. Whether to print registration success messages.
#'   Default: inherits from package option \code{getOption("SigBridgeRUtils.verbose")}.
#'
#' @return Invisibly returns \code{TRUE} on successful registration (via the underlying registrar).
#' @export
#' @keywords Registering
#'
#' @examples
#' \dontrun{
#' # Register a screening method for binary/survival phenotypes
#' Register(
#'   registry = "ScreenStrategy",
#'   Scissor = DoScissor,
#'   supported_phenotypes = c("binary", "survival")
#' )
#'
#' # Register a preprocessing step (e.g., Harmony integration)
#' Register(
#'   registry = "SCPreProcessStrategy",
#'   h = "Seurat::RunHarmony"
#' )
#'
#' # Register an annotation method
#' Register(
#'   registry = "AnnotationStrategy",
#'   my_annot = MyCustomAnnotator
#' )
#'
#' # Partial matching works
#' Register(registry = "anno", my_annot2 = AnotherAnnotator)
#' }
#'
#' @seealso
#'   \code{\link{RegisterScreenMethod}},
#'   \code{\link{RegisterSeuratMethod}},
#'   \code{\link{RegisterAnnoMethod}},
#'   \code{\link{SCPreProcessStrategy}}
Register <- function(
  ...,
  registry = c("ScreenStrategy", "SCPreProcessStrategy", "AnnotationStrategy"),
  verbose = getFuncOption("verbose")
) {
  #   dots <- rlang::list2(...)
  #   verbose <- dots$verbose %||% getFuncOption("verbose")
  if (is.character(registry)) {
    registry <- SigBridgeRUtils::MatchArg(
      registry,
      c("ScreenStrategy", "SCPreProcessStrategy", "AnnotationStrategy")
    )
    registry_obj <- get0(registry)
  } else {
    registry_obj <- registry
  }

  switch(
    registry,
    "ScreenStrategy" = RegisterScreenMethod(
      ...,
      registry = registry_obj,
      verbose = verbose
    ),
    "SCPreProcessStrategy" = RegisterSeuratMethod(
      ...,
      registry = registry_obj,
      verbose = verbose
    ),
    "AnnotationStrategy" = RegisterAnnoMethod(
      ...,
      registry = registry_obj,
      verbose = verbose
    )
  )
}
