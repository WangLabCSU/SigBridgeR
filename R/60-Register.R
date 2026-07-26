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
#'   \item \code{registry = "auto"} → \code{detect_registry()}
#'   \item \code{registry = "ScreenStrategy"} → \code{\link{RegisterScreenMethod}}
#'   \item \code{registry = "SCPreProcessStrategy"} → \code{\link{RegisterSeuratMethod}}
#'   \item \code{registry = "SCAnnotateStrategy"} → \code{\link{RegisterAnnoMethod}}
#' }
#'
#' @param ... Arguments passed to the underlying registrar. The exact requirements
#'   depend on the target \code{registry}:
#'   \describe{
#'     \item{\code{ScreenStrategy}}{Named functions with optional \code{supported_phenotypes},
#'       \code{parameter_mapper}, etc. (see \code{\link{RegisterScreenMethod}}).}
#'     \item{\code{SCPreProcessStrategy}}{Named functions or character specifications
#'       (e.g., \code{"h" = "Seurat::RunHarmony"}; see \code{\link{RegisterSeuratMethod}}).}
#'     \item{\code{SCAnnotateStrategy}}{Named annotation functions (see \code{\link{RegisterAnnoMethod}}).}
#'   }
#' @param registry Character. Target strategy environment for registration.
#'   Must be one of: \code{"auto"}, \code{"ScreenStrategy"}, \code{"SCPreProcessStrategy"}, or
#'   \code{"SCAnnotateStrategy"}. Partial matching is supported (e.g., \code{"screen"} → \code{"ScreenStrategy"}).
#' @param verbose Logical. Whether to print registration success messages.
#'   Default: inherits from package option \code{getOption("SigBridgeRUtils.verbose")}.
#'
#' @return Invisibly returns \code{TRUE} on successful registration (via the underlying registrar).
#' @export
#' @name Register
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
#'   registry = "SCAnnotateStrategy",
#'   my_annot = MyCustomAnnotator
#' )
#'
#' # auto detects the target registry
#' Register(my_annot2 = AnotherAnnotator)
#' }
#'
#' @seealso
#'   \code{\link{RegisterScreenMethod}},
#'   \code{\link{RegisterSeuratMethod}},
#'   \code{\link{RegisterAnnoMethod}},
#'   \code{\link{SCPreProcessStrategy}}
#'   \code{\link{SCAnnotateStrategy}}
#'   \code{\link{ScreenStrategy}}
Register <- new_generic(
  name = "Register",
  dispatch_args = c("method", "name")
)

method(generic = Register, class = class_any) <- function(method, name, ...) {
  cls_method <- class(method)
  expected_cls <- c("ScreenMethod", "AnnotationMethod")
  Abort(
    "Unsupported class: {.cls {cls_method}}",
    "Expected {.cls {expected_cls}}"
  )
}

#' @rdname Register
#' @export
method(generic = Register, class = ScreenMethod) <- function(
  func,
  name,
  verbose = getFuncOption("verbose")
) {
  registry <- unlist(registry)
  chk::chk_logical(verbose)
  if (!is.character(registry)) {
    chk::chk_environment(registry)
  }

  if (is.character(registry)) {
    registry <- SigBridgeRUtils::MatchArg(
      registry,
      c(
        "auto",
        "ScreenStrategy",
        "SCPreProcessStrategy",
        "SCAnnotateStrategy"
      )
    )

    if (identical(registry, "auto")) {
      dots <- rlang::list2(...)
      is_func <- purrr::map_lgl(dots, base::is.function)
      for (i in seq_len(sum(is_func))) {
        registry <- detect_registry(
          method_name = names(dots)[i],
          func = dots[[i]],
          dots = dots,
          verbose = verbose
        ) # character
        registry_obj <- utils::getFromNamespace(registry, "SigBridgeR")

        switch(
          registry,
          "ScreenStrategy" = RegisterScreenMethod(
            !!!dots[is_func][i],
            registry = registry_obj,
            verbose = verbose
          ),
          "SCPreProcessStrategy" = RegisterSeuratMethod(
            !!!dots[is_func][i],
            registry = registry_obj,
            verbose = verbose
          ),
          "SCAnnotateStrategy" = RegisterAnnoMethod(
            !!!dots[is_func][i],
            registry = registry_obj,
            verbose = verbose
          )
        )
      }

      return(invisible(TRUE))
    }

    # "ScreenStrategy",        "SCPreProcessStrategy",        "SCAnnotateStrategy"
    registry_obj <- utils::getFromNamespace(registry, "SigBridgeR")
  } else if (is.environment(registry)) {
    chk::chk_length(registry)
    registry_obj <- registry
  } else {
    cli::cli_abort(c(
      "x" = "Unknown registry provided",
      ">" = "Expected {.cls {c('character', 'environment','list')}}, \
      got {.cls {class(registry)}}"
    ))
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
    "SCAnnotateStrategy" = RegisterAnnoMethod(
      ...,
      registry = registry_obj,
      verbose = verbose
    )
  )

  # invisbly TRUE
}

#' @rdname Register
#' @export
method(generic = Register, class = AnnotationMethod) <- function(
  func,
  name,
  verbose = getFuncOption("verbose")
) {}

#' @rdname Register
#' @export
method(generic = Register, class = class_function) <- function(
  func,
  name,
  verbose = getFuncOption("verbose")
) {
  # Seurat
}
