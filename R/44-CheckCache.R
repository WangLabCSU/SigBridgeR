#' @title Check Cache Configuration Consistency
#'
#' @description
#' Validates whether the cached configuration matches the current parameters.
#' This function reads a cache metadata JSON file and compares it with the
#' provided screen method, phenotype class, label type, and parameters.
#'
#' @param path Character string specifying the path to the cache directory
#'   or directly to the `cache_config.json` file.
#' @param screen_method Character string indicating the screening method used.
#' @param phenotype_class Character vector specifying the phenotype class type.
#'   Must be one of `"binary"`, `"survival"`, or `"continuous"`.
#' @param label_type Character string specifying the label type. Defaults to
#'   the value of `screen_method`.
#' @param params List of parameters used for the screening method.
#' @param ... Additional arguments (must be empty; raises error if provided).
#'
#' @return Returns `invisible(TRUE)` if the cache configuration is consistent
#'   with the current parameters. Otherwise, aborts with an error message
#'   displaying the differences.
#'
#' @family cache_config
#' @export
#'
CheckCache <- new_generic("CheckCache", dispatch_args = "cache_config")

method(
  generic = CheckCache,
  class = ScreenMethodConfig
) <- CheckCache.ScreenMethodConfig <- function(
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


method(generic = CheckCache, class = ScreenMethodCache) <- function(
  cache_config,
  path,
  ...
) {
  CheckCache.ScreenMethodConfig(
    cache_config = cache_config@screen_method_config,
    path = path,
    ...
  )
}
