#' @title Resolve and Set Up Cache Directory
#'
#' @description
#' Resolves a user-specified path to the cache layer directory. In
#' `"save"` mode, the necessary directory structure is created
#' automatically. In `"load"` mode, an existing cache directory is
#' identified.
#'
#' @details
#' The function detects which cache layer the given path belongs to:
#' \itemize{
#'   \item **Cache layer**: The path directly contains a `cache_config.json`
#'     file. In save mode this triggers an error (recursive caching is not
#'     supported); in load mode the path is returned as-is.
#'   \item **Root layer**: The path's basename matches `{method_name}_res`
#'     (e.g., `"Scissor_res"`). In save mode a new cache subdirectory is
#'     created; in load mode [ChooseCache()] is called to select an existing
#'     cache.
#'   \item **Parent layer**: The path is a parent directory. A
#'     `{method_name}_res` subdirectory is created (save) or located (load)
#'     underneath it.
#' }
#'
#' In save mode, the cache directory is named
#' `{phenotype_class}_{timestamp}` (e.g.,
#' `"survival_20250730120000"`).
#'
#' @param path `character`. User-specified path to a cache, root, or parent
#'   directory.
#' @param cache_config A `ScreenMethodConfig` or `ScreenMethodCache` S7
#'   object. `method_name` and `phenotype_class` are extracted from this
#'   object to construct the cache directory path.
#' @param mode `character`. Either `"load"` or `"save"`. Default: `"load"`.
#' @param timestamp `character` or `NULL`. Optional timestamp string for the
#'   new cache directory in save mode. When `NULL`, defaults to
#'   `format(Sys.time(), "%Y%m%d%H%M%S")`.
#' @param ... Additional arguments (must be empty, checked by
#'   `rlang::check_dots_empty0()`).
#'
#' @returns A `character` string: the absolute path to the cache directory.
#'
#' @family cache_config
#' @export
#'
#' @examples
#' \dontrun{
#' config <- ScreenMethodConfig(
#'   method_name = "Scissor",
#'   param = list(alpha = 0.05)
#' )
#'
#' # Save: create a new cache under the root
#' CacheSetHere("Scissor_res", config, mode = "save")
#'
#' # Load: select an existing cache from the root
#' CacheSetHere("Scissor_res", config, mode = "load")
#'
#' # Load: point directly to a specific cache
#' CacheSetHere("Scissor_res/survival_207701011212", config, mode = "load")
#' }
CacheSetHere <- function(
  path,
  cache_config,
  mode = c("load", "save"),
  timestamp = NULL,
  ...
) {
  check_dots_empty()
  chk::chk_string(path)
  mode <- arg_match(mode)
  phenotype_class <- cache_config@phenotype_class
  screen_method <- cache_config@method_name

  root_dir_name <- paste0(screen_method, "_res")

  # -- resolve to absolute path -------------------------------------------
  path <- normalizePath(path, mustWork = FALSE)

  # -- detect layer -------------------------------------------------------
  has_cache_config <- file.exists(file.path(path, "cache_config.json"))

  layer <- if (has_cache_config) {
    "cache"
  } else if (basename(path) == root_dir_name) {
    "root"
  } else {
    "parent"
  }

  # -- dispatch -----------------------------------------------------------
  switch(
    mode,
    "save" = save_impl(path, layer, root_dir_name, phenotype_class, timestamp),
    "load" = load_impl(path, layer, root_dir_name)
  )
}

# ---------------------------------------------------------------------------
# Save implementation
# ---------------------------------------------------------------------------

save_impl <- function(
  path,
  layer,
  root_dir_name,
  phenotype_class,
  timestamp
) {
  if (layer == "cache") {
    Abort(
      "Recursive caching is not supported.",
      "{.path {path}} already contains {.file cache_config.json}.",
      "Specify a root-level or parent-level path."
    )
  }

  root_dir <- if (layer == "root") {
    path
  } else {
    # parent layer
    root_dir <- file.path(path, root_dir_name)
    if (!dir.exists(root_dir)) {
      dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)
    }
    root_dir
  }

  cache_name <- generate_cache_name(phenotype_class, timestamp)
  cache_dir <- file.path(root_dir, cache_name)

  if (dir.exists(cache_dir)) {
    Abort(
      "Cache directory already exists: {.path {cache_dir}}",
      "Choose a different path or delete the existing cache."
    )
  }

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  invisible(normalizePath(cache_dir, mustWork = FALSE))
}

# ---------------------------------------------------------------------------
# Load implementation
# ---------------------------------------------------------------------------

load_impl <- function(path, layer, root_dir_name) {
  if (layer == "root") {
    # root layer -- delegate to ChooseCache
    if (!dir.exists(path)) {
      Abort("Root directory not found: {.path {path}}")
    }
    ChooseCache(path)
  } else if (layer == "cache") {
    # cache layer -- return as-is
    if (!dir.exists(path)) {
      Abort("Cache directory not found: {.path {path}}")
    }
    path
  } else {
    # parent layer -- look for root underneath
    root_dir <- file.path(path, root_dir_name)
    if (!dir.exists(root_dir)) {
      Abort(
        "No cache root directory found under {.path {path}}.",
        "Expected a directory named {.val {root_dir_name}}."
      )
    }
    ChooseCache(root_dir)
  }
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

generate_cache_name <- function(phenotype_class, timestamp = NULL) {
  if (is.null(timestamp)) {
    timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
  }
  paste(phenotype_class, timestamp, sep = "_")
}
