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
CheckCache <- function(
  path,
  screen_method,
  phenotype_class = c("binary", "survival", "continuous"),
  label_type = screen_method,
  params,
  ...
) {
  rlang::check_dots_empty0()
  rlang::check_installed("jsonlite")
  if (!endsWith(path, "cache_config.json") && dir.exists(path)) {
    meta_json <- file.path(path, "cache_config.json")
  } else if (endsWith(path, "cache_config.json")) {
    meta_json <- path
  } else {
    cli::cli_abort(c(
      "x" = "Unsupported path: {.file {path}}",
      ">" = "`path` must be either a path to {.file cache_config.json} , or a folder containing it"
    ))
  }
  # * a list, NULL value became a NULL list
  cache_meta <- make_null_list_NULL(jsonlite::read_json(meta_json))

  inputs_meta <- list(
    config = list(
      screen_method = screen_method,
      phenotype_class = phenotype_class,
      label_type = label_type,
      params = params
    )
  )

  modified_config <- utils::modifyList(
    cache_meta,
    inputs_meta,
    keep.null = TRUE
  )
  check_res <- identical(modified_config, cache_meta)
  if (check_res) {
    return(invisible(TRUE))
  }

  diff_in_vec <- find_diff_in_2_lists(cache_meta, modified_config)

  cli::cli_alert_warning(cli::col_yellow(
    "Found slots different from the current parameters:"
  ))
  print(diff_in_vec)

  cli::cli_abort(c(
    "x" = "Cache config is not consistent with the current parameters"
  ))
}

#' @keywords internal
make_null_list_NULL <- function(x) {
  if (!is.list(x)) {
    return(x)
  }

  if (length(x) == 0L) {
    return(NULL)
  }

  lapply(x, make_null_list_NULL)
}

#' @keywords internal
find_diff_in_2_lists <- function(x, y) {
  if (identical(x, y)) {
    return(TRUE)
  }

  purrr::map2(
    x,
    y,
    ~ {
      if (is.list(.x) && is.list(.y)) {
        return(find_diff_in_2_lists(.x, .y))
      }
      if (identical(.x, .y)) {
        return(TRUE)
      }
      return(FALSE)
    }
  )
}
