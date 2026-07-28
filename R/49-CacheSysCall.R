#' @inheritParams CacheSetHere
#' @export
CacheSysCall <- function(
  mode = c("load", "save"),
  path,
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
    timestamp = timestamp,
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
    LoadCache(file = file.path(path, param_name))
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
  EnsureParentDir(path, .envir = caller_env())
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
