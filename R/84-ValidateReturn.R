ValidateReturn <- function(func, ...) {
  values <- .vr_collect_return_values(fn_body(func))

  if (length(values) == 0L) {
    values <- list(NULL)
  }

  is_list <- vapply(values, .vr_is_list_creation, logical(1L))
  has_slot <- vapply(values, .vr_list_has_scRNA_data, logical(1L))

  failed <- FALSE

  if (!all(is_list)) {
    validate_error("Return value is not a list")

    cli::cli_div(
      id = "missing_return",
      theme = list(div = list("margin-left" = 2L))
    )

    validate_explain_speaker(
      message = c(
        "recommended to be the first element of the return value",
        "Should be of class {.cls Seurat}"
      ),
      message_name = "scRNA_data",
      type = "error"
    )

    cli::cli_end(id = "missing_return")

    failed <- TRUE
  }

  if (any(is_list & !has_slot)) {
    validate_error("Return value does not have `scRNA_data` slot")

    cli::cli_div(
      id = "missing_return",
      theme = list(div = list("margin-left" = 2L))
    )

    validate_explain_speaker(
      message = c("Should be of class {.cls Seurat}"),
      message_name = "scRNA_data",
      type = "error"
    )

    cli::cli_end(id = "missing_return")

    failed <- TRUE
  }

  if (failed) {
    return(list(error = 1L))
  }

  validate_success("Return value is a list with `scRNA_data` slot")

  invisible()
}

# Collects every expression a function may return, following all `if` branches,
# the last expression of every block, and explicit `return()` calls.
# @param expr A function body or any sub-expression of it.
# @return A list of expressions (`NULL` stands for a missing return value).
# @keywords internal
.vr_collect_return_values <- function(expr) {
  if (!is_call(expr)) {
    return(list(expr))
  }

  fn_name <- safe_call_name(expr)

  if (is.null(fn_name)) {
    return(list(expr))
  }

  if (fn_name == "return") {
    return(list(if (length(expr) >= 2L) expr[[2L]] else NULL))
  }

  if (fn_name == "{") {
    exprs <- safe_call_args(expr)

    if (length(exprs) == 0L) {
      return(list(NULL))
    }

    early_exits <- .vr_explicit_returns(exprs[-length(exprs)])
    last_value <- .vr_collect_return_values(exprs[[length(exprs)]])

    return(c(early_exits, last_value))
  }

  if (fn_name == "if") {
    # * both branches are reachable, so both have to comply
    then_value <- if (length(expr) >= 3L) {
      .vr_collect_return_values(expr[[3L]])
    } else {
      list(NULL)
    }

    else_value <- if (length(expr) >= 4L) {
      .vr_collect_return_values(expr[[4L]])
    } else {
      list(NULL)
    }

    return(c(then_value, else_value))
  }

  list(expr)
}

# Collects the values of explicit `return()` calls found in `exprs`.
# @keywords internal
.vr_explicit_returns <- function(exprs) {
  values <- list()

  for (expr in exprs) {
    for (value in .vr_explicit_returns_in(expr)) {
      values[length(values) + 1L] <- list(value)
    }
  }

  values
}

# Recursively collects the values of explicit `return()` calls in `expr`.
# Nested function definitions are skipped, their `return()` belongs to them.
# @keywords internal
.vr_explicit_returns_in <- function(expr) {
  if (!is_call(expr)) {
    return(list())
  }

  fn_name <- safe_call_name(expr)

  if (is.null(fn_name) || fn_name %in% c("function", "~")) {
    return(list())
  }

  if (fn_name == "return") {
    return(list(if (length(expr) >= 2L) expr[[2L]] else NULL))
  }

  values <- list()

  for (arg in safe_call_args(expr)) {
    for (value in .vr_explicit_returns_in(arg)) {
      values[length(values) + 1L] <- list(value)
    }
  }

  values
}

# @keywords internal
.vr_is_list_creation <- function(expr) {
  is_call(expr) && identical(safe_call_name(expr), "list")
}

# @keywords internal
.vr_list_has_scRNA_data <- function(expr) {
  if (!.vr_is_list_creation(expr)) {
    return(FALSE)
  }

  args <- safe_call_args(expr)
  arg_names <- names(args)

  if (!is.null(arg_names) && "scRNA_data" %in% arg_names) {
    return(TRUE)
  }

  # * untagged arguments may still carry the tag inside a `=` call
  any(vapply(
    args,
    function(arg) {
      is_call(arg) &&
        identical(safe_call_name(arg), "=") &&
        length(arg) >= 3L &&
        identical(as.character(arg[[2L]]), "scRNA_data")
    },
    logical(1L)
  ))
}


#' @keywords internal
safe_call_name <- function(expr) {
  tryCatch(call_name(expr), error = function(e) NULL)
}

#' @keywords internal
safe_call_args <- function(expr) {
  tryCatch(call_args(expr), error = function(e) list())
}
