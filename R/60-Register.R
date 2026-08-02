#' @title Unified Registration Interface for Strategy Methods
#'
#' @description
#' Registers screening, preprocessing, or annotation methods into their
#' respective strategy registries. `Register()` accepts named arguments where
#' the name becomes the registry key and the value (an S7 object or function)
#' is dispatched to the appropriate `RegisterImpl` method.
#'
#' @details
#' `Register()` iterates over its named `...` arguments, calling
#' `RegisterImpl()` for each. `RegisterImpl` is an S7 generic that dispatches
#' based on the class of the value:
#' \itemize{
#'   \item [ScreenMethod] → registered into [ScreenStrategy].
#'   \item [AnnotationMethod] → registered into [SCAnnotateStrategy].
#'   \item `function` → validated as a Seurat/SeuratObject function and
#'     registered into [SCPreProcessStrategy].
#' }
#'
#' If the method name already exists in the target registry, the function
#' aborts unless `overwrite = TRUE`.
#'
#' @param ... Named arguments where each name is the registry key and each
#'   value is an S7 object ([ScreenMethod] or `AnnotationMethod`) or a
#'   function (for `SCPreProcessStrategy`).
#' @param overwrite `logical`. Whether to overwrite an existing entry with
#'   the same name. Default: `FALSE`.
#' @param verbose `logical`. Whether to print registration success messages.
#'   Default: `getFuncOption("verbose")`.
#'
#' @returns Invisibly `TRUE` on successful registration.
#'
#' @seealso
#'   [ScreenStrategy],
#'   [SCAnnotateStrategy],
#'   [SCPreProcessStrategy]
#'
#' @name Register
#' @export
#'
#' @examples
#' \dontrun{
#' # Register a ScreenMethod
#' method <- ScreenMethod(
#'   method_name = "Scissor",
#'   executor = DoScissor
#' )
#' Register(Scissor = method)
#'
#' # Register a Seurat preprocessing function
#' Register(h = Seurat::RunHarmony)
#'
#' # Register with overwrite
#' Register(Scissor = method, overwrite = TRUE)
#' }
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
  dispatch_args = "x"
)

method(generic = RegisterImpl, class_any) <- function(
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

method(generic = RegisterImpl, ScreenMethod) <- function(
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

method(generic = RegisterImpl, AnnotationMethod) <- function(
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

method(generic = RegisterImpl, class_function) <- function(
  x,
  name,
  overwrite = FALSE,
  verbose = getFuncOption("verbose")
) {
  # Seurat
  chk::chk_character(name, "name")
  if (nchar(name) != 1) {
    Abort(
      "Name (key) must be a single letter",
      "Current Name (key): {.val {name}}",
      "If you intend to register a screen method or annotation method,\
       please create a {.cls ScreenMethod} or {.cls AnnotationMethod} object"
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
