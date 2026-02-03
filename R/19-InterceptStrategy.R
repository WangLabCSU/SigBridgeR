#' @title Inspect Registered Strategy Environments
#'
#' @description
#' Lists and optionally inspects the contents of predefined strategy environments
#' (e.g., \code{ScreenStrategy}, \code{SCPreProcessStrategy}, \code{AnnotationStrategy}).
#' Returns a summary of registered methods or variables, with optional detailed
#' introspection of function signatures, body structure, and closure properties.
#'
#' This utility aids in debugging, validation, and exploration of dynamically
#' registered pipelines—particularly useful when developing or auditing extensible
#' analysis frameworks.
#'
#' @param strategy Character. Name of the strategy environment to inspect.
#'   Must be one of: \code{"ScreenStrategy"}, \code{"SCPreProcessStrategy"},
#'   or \code{"AnnotationStrategy"}. Partial matching is supported.
#' @param detailed Logical. If \code{TRUE}, performs deep inspection of function objects
#'   (e.g., argument count, presence of \code{...}, body length, closure environment).
#'   Default: \code{FALSE}.
#' @param verbose Logical. Whether to print a formatted table of results using \code{knitr::kable()}.
#'   Default: \code{TRUE}.
#' @param ... Additional arguments (currently unused, reserved for future extension).
#'
#' @return Invisibly returns a \code{data.frame}:
#'   \itemize{
#'     \item If \code{detailed = FALSE}: one row per entry in the strategy environment,
#'           with columns depending on content type (e.g., \code{method_name}, \code{func}, or \code{var}).
#'     \item If \code{detailed = TRUE}: augmented with introspection columns such as
#'           \code{func_n_args}, \code{func_has_dots}, \code{func_body_nlines},
#'           \code{func_is_closure}, etc.
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' # List all preprocessing steps
#' InterceptStrategy("SCPreProcessStrategy")
#'
#' # Get detailed function metadata for screen methods
#' details <- InterceptStrategy("ScreenStrategy", detailed = TRUE, verbose = FALSE)
#' head(details[, grep("executor_", names(details))])
#' }
InterceptStrategy <- function(
  strategy = c("ScreenStrategy", "SCPreProcessStrategy", "AnnotationStrategy"),
  detailed = FALSE,
  ...
) {
  strategy <- SigBridgeRUtils::MatchArg(
    strategy,
    c("ScreenStrategy", "SCPreProcessStrategy", "AnnotationStrategy"),
    NULL
  )

  vars <- names(get0(strategy))
  info <- dplyr::bind_rows(purrr::map(vars, function(var_name) {
    var <- get0(strategy)[[var_name]]
    if (length(var) > 1) {
      purrr::map(var, ~ if (is.function(.x)) list(.x) else .x)
    } else if (is.function(var)) {
      list(method_name = var_name, func = list(var))
    } else {
      list(method_name = var_name, var = var)
    }
  }))

  if (detailed) {
    detailed_info <- switch(
      strategy,
      "ScreenStrategy" = find_func_detail(info, "executor") %>%
        find_func_detail("mapper"),
      "SCPreProcessStrategy" = find_func_detail(info, "func"),
      "AnnotationStrategy" = find_func_detail(info, "func")
    )
  }

  if (detailed) {
    invisible(detailed_info)
  } else {
    invisible(info)
  }
}

#' @keywords internal
find_func_detail <- function(tbl, col_name) {
  if (!col_name %in% names(tbl)) {
    cli::cli_abort("Column {.val {col_name}} not found in tbl")
  }

  fn <- tbl[[col_name]]

  # 参数信息检查
  n_args <- purrr::map_int(
    fn,
    ~ {
      f <- formals(.x)
      if (is.null(f)) 0L else length(f)
    }
  )

  has_dots <- purrr::map_lgl(
    fn,
    ~ {
      f <- formals(.x)
      !is.null(f) && any(names(f) == "...")
    }
  )

  arg_names <- purrr::map_chr(
    fn,
    ~ {
      f <- formals(.x)
      if (is.null(f)) "" else paste(names(f), collapse = ", ")
    }
  )

  # 函数体结构检查
  body_length <- purrr::map_int(
    fn,
    ~ {
      b <- rlang::fn_body(.x)
      if (is.null(b) || is.primitive(.x)) 0L else length(b)
    }
  )

  body_nlines <- purrr::map_int(
    fn,
    ~ {
      b <- rlang::fn_body(.x)
      if (is.null(b) || is.primitive(.x)) 0L else length(deparse(b))
    }
  )

  # 环境与闭包特性检查
  is_closure <- purrr::map_lgl(
    fn,
    ~ {
      env <- environment(.x)
      !is.null(env) && !identical(env, emptyenv())
    }
  )

  env_parent <- purrr::map_chr(
    fn,
    ~ {
      env <- environment(.x)
      if (is.null(env)) {
        "NULL"
      } else if (identical(env, .GlobalEnv)) {
        "global"
      } else {
        "custom"
      }
    }
  )

  # 类型与来源检查
  is_primitive <- purrr::map_lgl(fn, base::is.primitive)

  is_s3_generic <- purrr::map_lgl(
    fn,
    ~ {
      b <- rlang::fn_body(.x)
      is.call(b) && as.character(b[[1]]) == "UseMethod"
    }
  )

  # 源码引用信息
  src_ref <- purrr::map_chr(
    fn,
    ~ {
      srcref <- attr(.x, "srcref")
      if (is.null(srcref) || is.primitive(.x)) {
        "unknown"
      } else {
        paste(srcref[1], collapse = ":")
      }
    }
  )

  # 构建带前缀的新列名
  prefix <- col_name
  new_cols <- list(
    n_args = n_args,
    has_dots = has_dots,
    arg_names = arg_names,
    body_length = body_length,
    body_nlines = body_nlines,
    is_closure = is_closure,
    env_parent = env_parent,
    is_primitive = is_primitive,
    is_s3_generic = is_s3_generic,
    src_ref = src_ref
  )
  # 获取原列位置并重组数据框（新列插入在原列之后）
  for (nm in names(new_cols)) {
    col_name_new <- paste0(prefix, "_", nm)
    tbl[[col_name_new]] <- new_cols[[nm]]
  }

  tbl
}
