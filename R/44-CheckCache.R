#' @title Check Cache Configuration Consistency
#'
#' @description
#' Validates whether the cached configuration matches the current parameters.
#' This function reads a cache metadata JSON file and compares it with the
#' provided screen method, phenotype class, label type, and parameters.
#'
#' # Methods
#' `CheckCache` is an S7 generic with methods available for the following
#' classes:
#'
#' `r doclisting::methods_list("CheckCache")`
#'
#' @details
#' The function compares each field in the current configuration against the
#' cached metadata. If any field differs (method name, method version,
#' phenotype class, label type, or parameters), the function aborts with a
#' formatted error message showing which parameter differs and both the
#' provided and actual values.
#'
#' @param cache_config A `ScreenMethodConfig` or `ScreenMethodCache` S7
#'   object containing the current configuration to validate.
#' @param ... Additional arguments to be passed to the method.
#'   - `path`: `character`. Path to the cache directory or directly to the
#'   `cache_config.json` file. Required for `ScreenMethodConfig`; derived
#'   from `cache_config@cache_config_path` for `ScreenMethodCache`.
#'
#' @returns Invisible `TRUE` if the cache configuration is consistent with
#'   the current parameters. Otherwise, aborts with an error message
#'   displaying the differences.
#'
#' @family cache_config
#' @name CheckCache
#' @export
#'
#' @examples
#' \dontrun{
#' # Check cache consistency for a ScreenMethodConfig
#' config <- ScreenMethodConfig(
#'   method_name = "Scissor",
#'   param = list(alpha = 0.05)
#' )
#' CheckCache(config, path = "cache/cache_config.json")
#'
#' # Check using a ScreenMethodCache object
#' cache <- ScreenMethodCache(
#'   cache_path = "cache/",
#'   cache_config_path = "cache/cache_config.json",
#'   cache_data = list(),
#'   screen_method_config = config
#' )
#' CheckCache(cache)
#' }
CheckCache <- new_generic("CheckCache", dispatch_args = "cache_config")

#' @rdname CheckCache
#' @export
method(CheckCache, class_any) <- function(
  cache_config,
  ...
) {
  cls_cache <- class(cache_config)
  expected_cls <- c("ScreenMethodConfig", "ScreenMethodCache")
  Abort(
    "cache_config must be a class of {.cls {expected_cls}}",
    "Current class is {.cls {cls_cache}}"
  )
}


CheckCache.ScreenMethodConfig <- function(
  cache_config,
  path,
  ...
) {
  check_installed("jsonlite")
  meta_json <- if (!endsWith(path, "cache_config.json") && dir.exists(path)) {
    file.path(path, "cache_config.json")
  } else if (endsWith(path, "cache_config.json")) {
    path
  } else {
    Abort(
      "Unsupported path: {.file {path}}",
      "`path` must be either a path to {.file cache_config.json} , or a folder containing it"
    )
  }
  # * a list, NULL value became a NULL list
  cache_meta <- make_null_list_NULL(jsonlite::read_json(meta_json)$config) # cpp fn
  # method_name
  # method_version
  # screen_method
  # phenotype_class
  # label_type
  # params = list(...)

  inputs_meta <- props(cache_config)

  diff <- find_diff_in_2_lists(
    inputs_meta, # user param
    cache_meta # actual cache config
  ) # cpp fn

  if (isTRUE(diff)) {
    return(invisible(TRUE))
  }

  arg_name <- diff$name
  user_val <- diff$x
  actual_val <- diff$y

  Abort(
    "Cache config is not consistent with the current parameters",
    "Parameter:{.field {arg_name}}, provided: {.val {user_val}}, actual: {.val {actual_val}}"
  )
}

#' @rdname CheckCache
#' @export
method(
  CheckCache,
  ScreenMethodConfig
) <- CheckCache.ScreenMethodConfig


#' @rdname CheckCache
#' @export
method(CheckCache, ScreenMethodCache) <- function(
  cache_config,
  ...
) {
  CheckCache.ScreenMethodConfig(
    cache_config = cache_config@screen_method_config,
    path = cache_config@cache_config_path
  )
}
