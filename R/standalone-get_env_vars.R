# ---
# repo: WangLabCSU/SigBridgeR
# file: standalone-get_env_vars.R
# last-updated: 2026-07-29
# license: https://unlicense.org
# dependencies: base
# standalone: yes
# ---
#
# This standalone file provides a lightweight helper to capture all variables
# (including hidden ones starting with `.`) from a given environment as a
# named list. It is self-contained (only uses `base` functions) and can be
# sourced independently or imported via `source()` / `box::use()` without
# pulling in additional package dependencies.
#
# Use cases include:
# - Debugging: snapshot the current local environment state.
# - Metaprogramming: capture all objects created inside a function's body.
# - Caching / reproducibility: serialise the environment to disk for later
#   inspection or replay.
#
# ## Functions
#
# - `get_env_vars(env)` – returns a named list of all bindings in `env`.
#
# ## Example
#
# ```r
# source("R/standalone-get_env_vars.R")
# f <- function() {
#   x <- 1
#   y <- "hello"
#   get_env_vars()
# }
# f()
# ```
#
# nocov start

#' Capture All Variables from an Environment
#'
#' @description
#' Captures every binding (including hidden variables whose names start
#' with `.`) from the specified environment and returns them as a **named
#' list**. The function only uses base-R facilities, making it suitable for
#' standalone inclusion in any script or package.
#'
#' @details
#' `get_env_vars()` calls `ls(envir = env, all.names = all.names)` to list all
#' names, then `mget(vars, envir = env, inherits = FALSE)` to retrieve the
#' corresponding values without walking up the parent chain.
#' `inherits = FALSE` ensures that only bindings directly defined in `env`
#' are returned, not those inherited from enclosing frames.
#'
#' The default `env = parent.frame()` captures the caller's environment,
#' which is the most common use-case inside a function body. If you need to
#' inspect another environment (e.g. `globalenv()`, a package namespace, or
#' a specific function's closure), pass it explicitly.
#'
#' @param env An environment whose bindings should be captured.
#'   Defaults to `parent.frame()`, i.e. the calling environment.
#'
#' @return A named list. Each element corresponds to one binding in `env`.
#'   Hidden names (starting with `.`) are included.
#'
#' @noRd
get_env_vars <- function(env = parent.frame(), all.names = FALSE) {
  vars <- ls(envir = env, all.names = all.names)
  setNames(mget(vars, envir = env, inherits = FALSE), vars)
}

# nocov end
