#' @title Write Cache Metadata
#'
#' @description
#' Writes metadata for SigBridgeR cache identification to a JSON file.
#' The metadata includes general system information, cache configuration,
#' and an optional user-provided description.
#'
#' # Methods
#' `WriteCacheMeta` is an S7 generic with methods available for the following
#' classes:
#'
#' `r doclisting::methods_list("WriteCacheMeta")`
#'
#' @details
#' This function creates a JSON metadata file containing:
#'
#' * **`general`**: File path, OS information, timestamp, R version, and
#'   SigBridgeR version
#' * **`config`**: Screen method, phenotype class, label type, and parameters
#' * **`description`**: Additional user-provided description
#'
#' The JSON file includes a header comment indicating it is required for
#' cache identification and should not be modified.
#'
#' Only for SigBridgeR version >= 3.7.0.
#'
#' @param cache_config A `ScreenMethodConfig` or `ScreenMethodCache` S7
#'   object containing the cache configuration to write.
#' @param ... Additional arguments passed to methods.
#'
#'   For `ScreenMethodConfig`:
#'   - `file`: `character(1)`. Path to the JSON metadata file to be written.
#'   - `additional_description`: `character` or `NULL`. An optional string
#'     describing the cache entry. Default: `NULL`.
#'   - `verbose`: `logical`. Whether to print progress messages. Default:
#'     `TRUE`.
#'
#'   For `ScreenMethodCache`:
#'   - `additional_description`: `character` or `NULL`. An optional string
#'     describing the cache entry. Default: `NULL`.
#'   - `verbose`: `logical`. Whether to print progress messages. Default:
#'     `TRUE`.
#'
#' @returns Invisible. Returns the metadata list that was written to the file.
#'
#' @family cache_config
#' @name WriteCacheMeta
#' @export
#'
#' @examples
#' \dontrun{
#' # Write metadata for a ScreenMethodConfig
#' config <- ScreenMethodConfig(
#'   method_name = "Scissor",
#'   param = list(alpha = 0.05)
#' )
#' WriteCacheMeta(config, file = "cache/cache_config.json")
#'
#' # With an additional description
#' WriteCacheMeta(config,
#'   file = "cache/cache_config.json",
#'   additional_description = "Scissor screening for survival phenotype"
#' )
#' }
WriteCacheMeta <- new_generic(
  name = "WriteCacheMeta",
  dispatch_args = "cache_config"
)

method(WriteCacheMeta, class_any) <- function(
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


WriteCacheMeta.ScreenMethodConfig <- function(
  cache_config,
  file,
  additional_description = NULL,
  verbose = TRUE,
  ...
) {
  check_dots_empty()
  check_installed("jsonlite")
  if (!is.null(additional_description)) {
    chk::chk_string(additional_description)
  }
  chk::chk_chr(file)

  cache_meta <- list(
    general = list2(
      file = file.path(getwd(), file),
      os = get_os_info(),
      time = get_time(),
      !!!get_r_info(),
      SigBridgeR_version = get_pkg_version()
    ),
    config = props(cache_config), # all properties as a list
    description = additional_description
  )

  if (verbose) {
    cli::cli_alert_info("Meta to: {.path {file}}")
  }

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

  json_content <- jsonlite::toJSON(
    cache_meta,
    auto_unbox = TRUE,
    pretty = TRUE
  )

  con <- file(file, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(json_content, con)

  if (verbose) {
    cli::cli_alert_success("Successfully wrote cache meta data.")
  }

  invisible(cache_meta)
}

method(WriteCacheMeta, ScreenMethodConfig) <- function(
  cache_config,
  ...,
  file,
  additional_description = NULL,
  verbose = TRUE
) {
  WriteCacheMeta.ScreenMethodConfig(
    cache_config = cache_config,
    file = file,
    additional_description = additional_description,
    verbose = verbose,
    ...
  )
}

method(WriteCacheMeta, ScreenMethodCache) <- function(
  cache_config,
  ...,
  additional_description = NULL,
  verbose = TRUE
) {
  WriteCacheMeta.ScreenMethodConfig(
    cache_config = cache_config@screen_method_config,
    file = cache_config@cache_config_path,
    additional_description = additional_description,
    verbose = verbose,
    ...
  )
}


get_os_info <- function() {
  sys_info <- Sys.info()
  os_info <- trimws(paste(sys_info["sysname"], sys_info["release"]))
  if (!nzchar(os_info)) {
    return(sys_info["sysname"])
  }
  os_info
}


get_r_info <- function() {
  list(
    r_version = R.version.string,
    arch = R.version$arch,
    platform = R.version$platform
  )
}


get_time <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}
