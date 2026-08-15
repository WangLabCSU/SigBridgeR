# ---
# repo: WangLabCSU/SigBridgeR
# file: standalone-get_var_value.R
# last-updated: 2026-07-18
# license: https://unlicense.org
# imports: [rlang]
# ---
#
# This file provides a tool to retrieve the value of a variable in a function.
#
# ## Changelog:
#
# 2026-07-18:
# * Changed internal label to noRd
#
# nocov start

# ---- Helper: initialise definition environment from formals ----------------
#' @noRd
.gvv_init_def_env <- function(fmls) {
  def_env <- new.env(parent = emptyenv())
  for (nm in names(fmls)) {
    if (nm == "...") {
      next
    } # ellipsis cannot be accessed via [[
    default <- fmls[[nm]]
    if (missing(default)) {
      next
    } # parameter without a default value
    if (!is.symbol(default) || nzchar(as.character(default))) {
      def_env[[nm]] <- default
    }
  }
  def_env
}

# ---- Helper: resolve expression by substituting known variables ------------
# Scans `expr` for symbols that exist in `def_env` and replaces them with
# their stored definition, walking the chain until literals or unresolvable
# symbols remain.
#' @noRd
.gvv_resolve <- function(expr, def_env, visited = character(), depth = 0L) {
  if (is.symbol(expr)) {
    name <- as.character(expr)
    if (name %in% visited) {
      rlang::abort(sprintf(
        "Circular dependency detected for variable: '%s'",
        name
      ))
    }
    if (exists(name, envir = def_env, inherits = FALSE)) {
      return(.gvv_resolve(
        def_env[[name]],
        def_env,
        c(visited, name),
        depth + 1L
      ))
    }
    return(expr)
  }
  if (is.call(expr)) {
    new_args <- vector("list", length(expr) - 1L)
    for (i in seq_along(expr)[-1L]) {
      new_args[[i - 1L]] <- .gvv_resolve(
        expr[[i]],
        def_env,
        visited,
        depth + 1L
      )
    }
    names(new_args) <- names(expr)[-1L]
    return(as.call(c(list(expr[[1L]]), new_args)))
  }
  expr # literal - return as-is
}

# ---- Helper: evaluate a fully-resolved expression safely ------------------
# Uses the evaluation environment stored inside `def_env` (key `.gvv_env`).
# Falls back to `baseenv()` when not set.
#' @noRd
.gvv_eval_safe <- function(expr, def_env) {
  eval_env <- def_env[[".gvv_env"]] %||% baseenv()
  tryCatch(eval(expr, envir = eval_env), error = function(e) NULL)
}

# ---- Helper: extract assignments from a `{` block -------------------------
# Processes statements sequentially; stops at `return()` / `stop()` or
# at `break` (signalled via the `.loop_break` flag in def_env).
#' @noRd
.gvv_extract_block <- function(expr, def_env) {
  for (i in seq_along(expr)[-1L]) {
    child <- expr[[i]]
    .gvv_extract(child, def_env)
    # Stop processing siblings if a terminating call or break was hit.
    # Two cases: plain symbol (return/stop/abort/cli_abort)
    # or namespace-qualified (rlang::abort / cli::cli_abort).
    if (is.call(child) && !is.call(child[[1L]])) {
      child_op <- as.character(child[[1L]])
      if (child_op %in% c("return", "stop", "abort", "cli_abort")) break
    } else if (is.call(child) && is.call(child[[1L]])) {
      ns_call <- child[[1L]]
      if (
        length(ns_call) >= 3L &&
          as.character(ns_call[[1L]]) == "::" &&
          paste0(
            as.character(ns_call[[2L]]),
            "::",
            as.character(ns_call[[3L]])
          ) %in%
            c("rlang::abort", "cli::cli_abort")
      ) {
        break
      }
    }
    if (isTRUE(def_env$.loop_break)) break
  }
}

# ---- Loop: `for (var in seq) body` ----------------------------------------
#' @noRd
.gvv_extract_for <- function(expr, def_env) {
  loop_var <- as.character(expr[[2L]])
  seq_expr <- .gvv_resolve(expr[[3L]], def_env)
  body_expr <- expr[[4L]]

  seq_value <- .gvv_eval_safe(seq_expr, def_env)

  if (
    !is.null(seq_value) &&
      is.atomic(seq_value) &&
      length(seq_value) > 0L &&
      length(seq_value) < 1e5L
  ) {
    for (iter_val in seq_value) {
      def_env[[loop_var]] <- iter_val
      .gvv_extract(body_expr, def_env)
    }
  } else {
    # Fallback: process the body once
    .gvv_extract(body_expr, def_env)
  }
}

# ---- Loop: `while (cond) body` --------------------------------------------
#' @noRd
.gvv_extract_while <- function(expr, def_env) {
  cond_expr <- expr[[2L]]
  body_expr <- expr[[3L]]

  for (iter in seq_len(10000L)) {
    cond_resolved <- .gvv_resolve(cond_expr, def_env)
    cond_value <- .gvv_eval_safe(cond_resolved, def_env)
    if (!isTRUE(cond_value)) {
      break
    }

    .gvv_extract(body_expr, def_env)
    if (isTRUE(def_env$.loop_break)) {
      def_env$.loop_break <- NULL
      return(invisible())
    }
    if (iter == 10000L) {
      rlang::warn("Reached iteration limit (10000) for while loop")
    }
  }
}

# ---- Loop: `repeat body` --------------------------------------------------
#' @noRd
.gvv_extract_repeat <- function(expr, def_env) {
  body_expr <- expr[[2L]]

  for (iter in seq_len(10000L)) {
    .gvv_extract(body_expr, def_env)
    if (isTRUE(def_env$.loop_break)) {
      def_env$.loop_break <- NULL
      return(invisible())
    }
    if (iter == 10000L) {
      rlang::warn(
        "Reached iteration limit (10000) for repeat loop, result may be unpredictable"
      )
    }
  }
}

# ---- Branch: `if (cond) then else else` -----------------------------------
# Evaluates the condition statically when possible; otherwise falls back
# to processing both branches.
#' @noRd
.gvv_extract_if <- function(expr, def_env) {
  cond_resolved <- .gvv_resolve(expr[[2L]], def_env)
  cond_value <- .gvv_eval_safe(cond_resolved, def_env)

  if (isTRUE(cond_value)) {
    .gvv_extract(expr[[3L]], def_env)
  } else if (identical(cond_value, FALSE)) {
    # FALSE with no else -> do nothing
    if (length(expr) >= 4L) .gvv_extract(expr[[4L]], def_env)
  } else {
    # Cannot determine branch statically - process both (conservative)
    .gvv_extract(expr[[3L]], def_env)
    if (length(expr) >= 4L) .gvv_extract(expr[[4L]], def_env)
  }
}

# ---- Main extractor dispatcher --------------------------------------------
# Walks the expression tree and records every assignment in `def_env`.
# Dispatches to specialised handlers for control-flow constructs so that
# iterative / conditional assignments are applied correctly.
#' @noRd
.gvv_extract <- function(expr, def_env) {
  if (!is.recursive(expr) || !is.call(expr)) {
    return(invisible())
  }

  # Function position is itself a call - namespace-qualified (pkg::fun)
  # or anonymous ((function(x) ...)(y)).
  if (is.call(expr[[1L]])) {
    ns_call <- expr[[1L]]
    # Check for terminating calls via namespace: rlang::abort, cli::cli_abort
    if (
      length(ns_call) >= 3L &&
        as.character(ns_call[[1L]]) == "::"
    ) {
      pkg <- as.character(ns_call[[2L]])
      fun <- as.character(ns_call[[3L]])
      if (paste0(pkg, "::", fun) %in% c("rlang::abort", "cli::cli_abort")) {
        return(invisible())
      }
    }
    # Other namespace calls - no assignments inside
    return(invisible())
  }

  op <- as.character(expr[[1L]])

  # --- Control-transfer statements ---
  if (op %in% c("return", "stop", "abort", "cli_abort")) {
    return(invisible())
  }
  if (op == "break") {
    def_env$.loop_break <- TRUE
    return(invisible())
  }
  if (op == "next") {
    return(invisible())
  } # skip rest of this iteration

  # --- Assignment operators ---
  if (op %in% c("<-", "=", "<<-")) {
    lhs_expr <- expr[[2L]]

    # Subset-assignments like x[i, j] <- value are represented as
    #   (`<-`, (`[`, x, i, j), value)
    # We transform them into a call to the replacement function `[<-`:
    #   x <- `[<-`(x, i, j, value = value)
    # and store the whole replacement expression under the base variable.
    if (
      is.call(lhs_expr) &&
        as.character(lhs_expr[[1L]]) %in% c("[", "[[", "$")
    ) {
      lhs <- as.character(lhs_expr[[2L]])
      # Collect subscripts (everything after the base object)
      subscripts <- as.list(lhs_expr)[-(1L:2L)]
      # Build: `[<-`(x, ..., value = RHS)
      replacement <- as.call(c(
        list(as.symbol(paste0(as.character(lhs_expr[[1L]]), "<-"))),
        lhs_expr[[2L]], # base object (x)
        subscripts, # indices (i, j, ...)
        expr[[3L]] # RHS value
      ))
      def_env[[lhs]] <- .gvv_resolve(replacement, def_env)
    } else {
      lhs <- as.character(lhs_expr)
      def_env[[lhs]] <- .gvv_resolve(expr[[3L]], def_env)
    }
    # fall through so that nested assignments in the RHS are also caught
  }

  # --- `{` block ---
  if (op == "{") {
    .gvv_extract_block(expr, def_env)
    return(invisible())
  }

  # --- Control-flow constructs ---
  switch(
    op,
    "for" = .gvv_extract_for(expr, def_env),
    "while" = .gvv_extract_while(expr, def_env),
    "repeat" = .gvv_extract_repeat(expr, def_env),
    "if" = .gvv_extract_if(expr, def_env),
    "::" = NULL, # namespace access - keep intact
    ":::" = NULL,
    {
      # Default: recursively traverse child expressions
      for (i in seq_along(expr)[-1L]) {
        .gvv_extract(expr[[i]], def_env)
      }
    }
  )
}

# ---- Public API -----------------------------------------------------------

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
#' @param .env Evaluation environment used for resolving function calls
#'   (e.g. `runif` from the **stats** package). Defaults to the caller's
#'   environment (`rlang::current_env()`), which is typically `.GlobalEnv`
#'   and has access to all attached packages.
#'
#' @return The computed value of the variable.
#'
#' @examples
#' f <- function(save_path = "./analysis") {
#'   save_path_new <- file.path(save_path, "res")
#'   return(save_path_new)
#' }
#' get_var_value("save_path_new", f)   # returns "./analysis/res"
#'
#' g <- function(a = 1, b = 2){
#'   c <- a * 2 + b * 3
#'   d = c^2
#'   e <<- d - 1
#'   e
#' }
#' get_var_value("e", g)   # returns 63
#'
#' h <- function(a = "A", ...){
#'   a <- 1
#'   return(a)
#'   a <- 2
#'   a
#' }
#' get_var_value("a",h)    # returns 1 (dead code after return is ignored)
#'
#' i <- function(x = 2) {
#'   for (k in 1:3) x <- x * 2
#'   x
#' }
#' get_var_value("x", i)   # returns 16
#'
#' j <- function(x = 2) {
#'   while (x < 10) x <- x * 2
#'   x
#' }
#' get_var_value("x", j)   # returns 16
#'
#' @noRd
get_var_value <- function(var_name, func, .env = rlang::current_env()) {
  # 1. Build definition environment from function formals
  def_env <- .gvv_init_def_env(formals(func))
  def_env[[".gvv_env"]] <- .env # store for .gvv_eval_safe

  # 2. Walk the function body and record all assignments
  .gvv_extract(rlang::fn_body(func), def_env)

  # 3. Verify the target variable exists
  if (!exists(var_name, envir = def_env, inherits = FALSE)) {
    rlang::abort(c(
      sprintf("Variable '%s' not found in function definition.", var_name),
      ">" = "It is neither a parameter with default nor assigned in the function body."
    ))
  }

  # 4. Resolve the target's expression (substitute known variables)
  resolved <- .gvv_resolve(def_env[[var_name]], def_env)

  # 5. Evaluate and return
  eval(resolved, envir = .env)
}
