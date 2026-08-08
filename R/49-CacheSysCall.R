#' @title System-Level Cache Load and Save
#'
#' @description
#' High-level wrapper for loading and saving cached screening results.
#' Combines cache path resolution (via [CacheSetHere()]), configuration
#' validation (via [CheckCache()]), and data I/O (via [LoadCache()] or
#' [WriteCache()]) into a single call.
#'
#' @details
#' `CacheSysCall` is the primary entry point for cache operations in the
#' SigBridgeR pipeline. It dispatches to either `CacheSysCall.load()` or
#' `CacheSysCall.save()` based on the `mode` argument.
#'
#' ## Load mode
#'
#' In `"load"` mode, the function:
#'
#' 1. Resolves the cache path via [CacheSetHere()].
#' 2. Validates cache consistency via [CheckCache()].
#' 3. Loads each cached data object via [LoadCache()] and returns them
#'    as a named list.
#'
#' ## Save mode
#'
#' In `"save"` mode, the function:
#'
#' 1. Resolves and creates the cache directory via [CacheSetHere()].
#' 2. Writes the cache data and metadata via [WriteCache()].
#' 3. Returns the cache data invisibly.
#'
#' @param mode `character`. Either `"load"` or `"save"`.
#' @param path `character` or `NULL`. User-specified path. In save mode,
#'   defaults to `cache@cache_path`. In load mode, defaults to `NULL` (user
#'   will be prompted to choose).
#' @param cache A `ScreenMethodConfig` or `ScreenMethodCache` S7 object.
#'   For save mode, a `ScreenMethodCache` is required (contains the data to
#'   write). For load mode, either class is accepted.
#' @param verbose `logical`. Whether to print progress messages.
#'   Default: `TRUE`.
#' @param timestamp `character` or `NULL`. Optional timestamp string for the
#'   new cache directory in save mode. Passed to [CacheSetHere()].
#' @param ... Additional arguments passed to [WriteCache()] in save mode
#'   (e.g., `format`, `max_rows`, `max_cols`).
#'
#' @returns
#' In `"load"` mode, a named list of loaded cache data objects. In
#' `"save"` mode, the cache data invisibly.
#'
#' @family cache_config
#' @export
#'
#' @examples
#' \dontrun{
#' # Load cached results
#' cache <- ScreenMethodCache(
#'   cache_path = "Scissor_res/",
#'   cache_config_path = "Scissor_res/cache_config.json",
#'   cache_data = list(),
#'   screen_method_config = ScreenMethodConfig(
#'     screen_method = "Scissor",
#'     param = list(alpha = 0.05)
#'   )
#' )
#' result <- CacheSysCall("load", cache = cache)
#'
#' # Save results to cache
#' cache <- ScreenMethodCache(
#'   cache_path = "Scissor_res/",
#'   cache_config_path = "Scissor_res/cache_config.json",
#'   cache_data = list(result = my_matrix),
#'   screen_method_config = ScreenMethodConfig(
#'     screen_method = "Scissor",
#'     param = list(alpha = 0.05)
#'   )
#' )
#' CacheSysCall("save", cache = cache)
#' }
CacheSysCall <- function(
  mode = c("load", "save"),
  path = if (mode == "save") cache@cache_path else NULL,
  cache,
  verbose = TRUE,
  timestamp = NULL,
  ...
) {
  mode <- arg_match(mode)
  cache_config <- if (S7_inherits(cache, "ScreenMethodConfig")) {
    if (mode == "save") {
      Abort(
        "cache is a {.cls ScreenMethodConfig}, but mode is {.val save}.",
        "Try passing a {.cls ScreenMethodCache} instead."
      )
    }
    cache
  } else if (S7_inherits(cache, "ScreenMethodCache")) {
    cache@screen_method_config
  } else {
    Abort(
      "cache is not a {.cls ScreenMethodConfig} or {.cls ScreenMethodCache}"
    )
  }

  path <- CacheSetHere(
    path = path,
    cache_config = cache_config,
    mode = mode,
    timestamp = timestamp
  )

  dots <- list2(...)

  switch(
    mode,
    # * Return a list containing the loaded data
    load = CacheSysCall.load(
      path = path,
      cache_config = cache_config,
      verbose = verbose
    ),
    # * Return a list containing the saved data invisibly
    save = exec(
      CacheSysCall.save,
      path = path,
      cache = cache,
      verbose = verbose,
      !!!dots
    )
  )
}

CacheSysCall.load <- function(
  path,
  cache_config,
  verbose = TRUE
) {
  CheckCache(cache_config = cache_config, path = path)
  param_names <- names(cache_config)

  if (length(param_names) == 0L) {
    Abort("param_names is empty, nothing to load")
  }

  load_fn <- purrr::in_parallel(\(param_name) {
    file <- file.path(path, param_name)
    if (!file.exists(file)) {
      return(NULL)
    }
    LoadCache(file = file)
  })
  purrr::map(
    param_names,
    load_fn,
    .progress = if (verbose) {
      "Loading cache"
    } else {
      FALSE
    }
  )
  # A LIST
}

CacheSysCall.save <- function(path, cache, verbose = TRUE, ...) {
  EnsureParentDir(path)
  WriteCache(
    cache,
    # format = c("auto", "qs2", "qdata", "rdata", "csv"),
    # max_rows = 1000L,
    # max_cols = 20L,
    # additional_description = NULL,
    verbose = verbose,
    ...
  )
}
