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
Register <- function(
  ...,
  overwrite = FALSE,
  verbose = getFuncOption("verbose")
) {
  dots <- list(...)
  purrr::iwalk(dots, \(x, name) {
    RegisterImpl(x = x, name = name, overwrite = overwrite, verbose = verbose)
  })
}

RegisterImpl <- new_generic(
  name = "RegisterImpl",
  dispatch_args = c("x", "name")
)


method(generic = RegisterImpl, class = class_any) <- function(
  x,
  name,
  ...
) {
  cls_method <- class(x)
  expected_cls <- c("ScreenMethod", "AnnotationMethod", "function")
  Abort(
    "Unsupported class: {.cls {cls_method}}",
    "Expected {.cls {expected_cls}}"
  )
}

method(generic = RegisterImpl, class = ScreenMethod) <- function(
  x,
  name = NULL,
  overwrite = FALSE,
  verbose = getFuncOption("verbose")
) {
  # x is ScreenMethod
  # method_name = property_chr,
  # method_version = property_chr,
  # executor = property_fn
  # phenotype_class = property_phenotype_class,
  # mapper = property_mapper_fn

  if (is.null(name)) {
    if (verbose) {
      cli::cli_alert_info(
        "Name (key) not provided, using {.field @method_name}: {.val {x@method_name}}"
      )
    }
  } else {
    x@method_name <- name
  }

  if (x@method_name %chin% names(ScreenStrategy) && !overwrite) {
    Abort(
      "Method {.field {x@method_name}} already exists",
      "Please use `overwrite = {.val TRUE}` to overwrite"
    )
  }

  ScreenStrategy[[name]] <- x

  if (verbose) {
    cli::cli_alert_success(
      "Registered {.field {x@method_name}} to {.cls ScreenStrategy}"
    )
  }

  invisible(TRUE)
}

method(generic = RegisterImpl, class = AnnotationMethod) <- function(
  x,
  name,
  overwrite = FALSE,
  verbose = getFuncOption("verbose")
) {
  if (is.null(name)) {
    if (verbose) {
      cli::cli_alert_info(
        "Name (key) not provided, using {.field @method_name}: {.val {x@method_name}}"
      )
    }
  } else {
    x@method_name <- name
  }

  if (x@method_name %chin% names(SCAnnotateStrategy) && !overwrite) {
    Abort(
      "Method {.field {x@method_name}} already exists",
      "Please use `overwrite = {.val TRUE}` to overwrite"
    )
  }
  SCAnnotateStrategy[[name]] <- x
  if (verbose) {
    cli::cli_alert_success(
      "Registered {.field {x@method_name}} to {.cls SCAnnotateStrategy}"
    )
  }

  invisible(TRUE)
}

method(generic = RegisterImpl, class = class_function) <- function(
  x,
  name,
  overwrite = FALSE,
  verbose = getFuncOption("verbose")
) {
  # Seurat
  chk::chk_character(name, "name")
  if (nchar(name) != 1) {
    Abort(
      "Name (key) must be a single character",
      "Current value: {.val {name}}"
    )
  }

  if (!is_func_from_pkg(x, "Seurat") && !is_func_from_pkg(x, "SeuratObject")) {
    Abort(
      "Function is not a Seurat function",
      "Please check typo",
      "Package version:\
       {.pkg Seurat}: {r_pkg_version('Seurat')},\
       {.pkg SeuratObject}: {r_pkg_version('SeuratObject')}"
    )
  }

  if (name %chin% names(SCPreProcessStrategy) && !overwrite) {
    Abort(
      "Method {.field {name}} already exists",
      "Please use `overwrite = {.val TRUE}` to overwrite"
    )
  }

  SCPreProcessStrategy[[name]] <- x
  if (verbose) {
    cli::cli_alert_success(
      "Registered {.field {name}} to {.cls SCPreProcessStrategy}"
    )
  }

  invisible(TRUE)
}

is_func_from_pkg <- function(func, pkg_name) {
  chk::chk_function(func)
  chk::chk_character(pkg_name)

  ns_name <- paste0("namespace:", pkg_name)

  env <- environment(func)

  # 快速路径：函数环境就是目标包 namespace
  if (!is.null(env) && identical(environmentName(env), ns_name)) {
    return(TRUE)
  }

  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    stop("pkg not installed")
  }

  ns <- asNamespace(pkg_name)

  # 若函数环境就是 namespace，也认为来自该包
  if (!is.null(env) && identical(env, ns)) {
    return(TRUE)
  }

  # 慢路径：在 namespace 中查找是否存在同一函数对象
  objs <- ls(ns, all.names = TRUE)

  for (nm in objs) {
    obj <- get(nm, envir = ns, inherits = FALSE)
    if (is.function(obj) && identical(obj, func)) {
      return(TRUE)
    }
  }

  FALSE
}
