#' @title Resolve and Set Up Cache Directory
#'
#' @description
#' Resolves a user-specified path to the cache layer directory. In
#' \code{"save"} mode, the necessary directory structure is created
#' automatically. In \code{"load"} mode, an existing cache directory is
#' identified.
#'
#' @param path Character string. User-specified path.
#' @param screen_method Character string. Screening method name. Must match
#'   a key in \code{\link{ScreenStrategy}}.
#' @param phenotype_class Character string. Phenotype class. One of
#'   \code{"binary"}, \code{"survival"}, or \code{"continuous"}.
#' @param mode Character string. Either \code{"load"} (default) or
#'   \code{"save"}.
#' @param timestamp Character string. Optional timestamp string for the new
#'   cache directory in save mode. Defaults to
#'   \code{format(Sys.time(), "\%Y_\%m_\%d_\%H\%M")}.
#' @param ... Additional arguments (must be empty).
#'
#' @return Character string. The absolute path to the cache layer directory.
#'
#' @family cache_config
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Save: create a new cache under the root
#' CacheSetHere("Scissor_res", "Scissor", "survival", mode = "save")
#'
#' # Load: select an existing cache from the root
#' CacheSetHere("Scissor_res", "Scissor", "survival", mode = "load")
#'
#' # Load: point directly to a specific cache
#' CacheSetHere("Scissor_res/survival_202512011212", "Scissor", "survival", mode = "load")
#' }
CacheSetHere <- function(
  path,
  screen_method,
  phenotype_class = c("binary", "survival", "continuous"),
  mode = c("load", "save"),
  timestamp = NULL,
  ...
) {
  rlang::check_dots_empty0()
  chk::chk_string(path)
  mode <- rlang::arg_match(mode)
  phenotype_class <- rlang::arg_match(phenotype_class)
  chk::chk_length(phenotype_class)
  screen_method <- rlang::arg_match(screen_method, names(ScreenStrategy))
  chk::chk_length(screen_method)

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

#' @keywords internal
save_impl <- function(
  path,
  layer,
  root_dir_name,
  phenotype_class,
  timestamp
) {
  if (layer == "cache") {
    cli::cli_abort(c(
      "x" = "Recursive caching is not supported.",
      "x" = "{.path {path}} already contains {.file cache_config.json}.",
      ">" = "Specify a root-level or parent-level path."
    ))
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
    cli::cli_abort(c(
      "x" = "Cache directory already exists: {.path {cache_dir}}",
      ">" = "Choose a different path or delete the existing cache."
    ))
  }

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  invisible(normalizePath(cache_dir, mustWork = FALSE))
}

# ---------------------------------------------------------------------------
# Load implementation
# ---------------------------------------------------------------------------

#' @keywords internal
load_impl <- function(path, layer, root_dir_name) {
  if (layer == "root") {
    # root layer -- delegate to ChooseCache
    if (!dir.exists(path)) {
      cli::cli_abort(c("x" = "Root directory not found: {.path {path}}"))
    }
    ChooseCache(path)
  } else if (layer == "cache") {
    # cache layer -- return as-is
    if (!dir.exists(path)) {
      cli::cli_abort(c("x" = "Cache directory not found: {.path {path}}"))
    }
    path
  } else {
    # parent layer -- look for root underneath
    root_dir <- file.path(path, root_dir_name)
    if (!dir.exists(root_dir)) {
      cli::cli_abort(c(
        "x" = "No cache root directory found under {.path {path}}.",
        ">" = "Expected a directory named {.val {root_dir_name}}."
      ))
    }
    ChooseCache(root_dir)
  }
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' @keywords internal
generate_cache_name <- function(phenotype_class, timestamp = NULL) {
  if (is.null(timestamp)) {
    timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
  }
  paste(phenotype_class, timestamp, sep = "_")
}
