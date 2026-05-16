# ---
# repo: WangLabCSU/SigBridgeR
# file: standalone-get_var_value.R
# last-updated: 2026-05-15
# license: https://unlicense.org
# imports: [rlang]
# ---

#' Trace and compute the value of a variable defined inside a function
#'
#' @description
#' Recursively traces variable assignments (both from function parameters and
#' function body) to reconstruct and evaluate the expression that produces
#' the variable's value. The trace follows the assignment chain until it
#' reaches literal values or unresolvable symbols.
#'
#' @param var_name Character string, the name of the target variable.
#' @param func An R function inside which the variable is defined.
#'
#' @return The computed value of the variable.
#'
#' @examples
#' f <- function(save_path = "./analysis") {
#'   save_path_new <- file.path(save_path, "res")
#'   return(save_path_new)
#' }
#' get_var_value("save_path_new", f) # returns "./analysis/res"
#'
#' g <- function(a = 1, b = 2){
#'   c <- a * 2 + b * 3
#'   d = c^2
#'   e <<- d - 1
#'   e
#' }
#' get_var_value("e", g) # returns 63
#'
#' h <- function(a = "A", ...){
#'   a <- 1
#'   return(a)
#'   a <- 2
#'   a
#' }
#' get_var_value("a",h) # returns 1
#'
#' j <- function(a = 1){
#'    print(a)
#'    a <- 2
#'    print(a)
#'    a
#' }
#'
#' @export
get_var_value <- function(var_name, func) {
  # ---- 1. Build variable-definition map ----
  fmls <- formals(func)
  body <- rlang::fn_body(func)

  # `var_defs` maps variable names to their defining R expressions
  var_defs <- list()

  # 1a. Extract parameter default values
  for (nm in names(fmls)) {
    if (nm == "...") {
      next
    } # skip ellipsis — cannot be accessed via [[
    default <- fmls[[nm]]

    if (missing(default)) {
      next
    }

    # Parameters without a default are empty symbols (name = "")
    if (!is.symbol(default) || nzchar(as.character(default))) {
      var_defs[[nm]] <- default
    }
  }

  # ---- 2. Resolve a (sub-)expression by substituting known variables ----
  #     (defined before extraction so assignment extraction can pre-resolve RHS)
  resolve_expr <- function(expr, visited = character()) {
    if (is.symbol(expr)) {
      name <- as.character(expr)
      if (name %in% names(var_defs)) {
        if (name %in% visited) {
          stop("Circular dependency detected for variable: ", name)
        }
        return(resolve_expr(var_defs[[name]], c(visited, name)))
      }
      return(expr) # keep as-is — final eval() will look it up
    }
    if (is.call(expr)) {
      new_args <- vector("list", length(expr) - 1L)
      for (i in seq_along(expr)[-1L]) {
        new_args[[i - 1L]] <- resolve_expr(expr[[i]], visited)
      }
      names(new_args) <- names(expr)[-1L]
      return(as.call(c(list(expr[[1L]]), new_args)))
    }
    expr # literal
  }

  # 1b. Extract assignments from the function body.
  #     **Pre-resolve** the RHS of each assignment immediately, so that
  #     self-referential assignments (e.g. `x <- x + 1` in a function whose
  #     parameter `x` has a default) use the OLD definition and avoid
  #     circular-dependency errors.
  #     Also stop at `return()` / `stop()` to ignore dead code.
  extract_assignments <- function(expr) {
    if (!is.recursive(expr)) {
      return(invisible())
    }
    if (is.call(expr)) {
      op <- as.character(expr[[1L]])

      if (op %in% c("return", "stop")) {
        return(invisible())
      }

      if (op %in% c("<-", "=", "<<-")) {
        lhs <- as.character(expr[[2L]])
        rhs <- expr[[3L]]
        # Resolve BEFORE storing so `x <- x + 1` becomes `x <- <old> + 1`
        var_defs[[lhs]] <<- resolve_expr(rhs)
      }

      if (op == "{") {
        for (i in seq_along(expr)[-1L]) {
          child <- expr[[i]]
          if (
            is.call(child) &&
              as.character(child[[1L]]) %in% c("return", "stop")
          ) {
            break
          }
          extract_assignments(child)
        }
      } else {
        for (i in seq_along(expr)[-1L]) {
          extract_assignments(expr[[i]])
        }
      }
    }
  }

  extract_assignments(body)

  if (!var_name %in% names(var_defs)) {
    stop(
      "Variable '",
      var_name,
      "' not found in function definition.\n",
      "It is neither a parameter with default nor assigned in the function body."
    )
  }

  # ---- 3. Resolve the target variable's expression ----
  #     For most cases the expression is already fully resolved (pre-resolved
  #     during extraction). This second pass handles cases where the target
  #     was not re-assigned in the body and still contains raw symbols.
  resolved <- resolve_expr(var_defs[[var_name]])

  # ---- 3. Evaluate the fully resolved expression ----
  eval(resolved, envir = baseenv())
}
