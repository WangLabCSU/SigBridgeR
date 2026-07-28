#' @title Write Cache Metadata
#'
#' @description Writes metadata for SigBridgeR cache identification to a JSON file.
#'
#' @param file Character string. Path to the JSON metadata file to be written.
#' @param ... Additional arguments (must be empty, checked by `rlang::check_dots_empty0()`).
#'
#' @return Invisible. Returns the metadata list that was written to the file.
#'
#' @details
#' This function creates a JSON metadata file containing:
#' \describe{
#'   \item{general}{File path, OS information, timestamp, R version, and SigBridgeR version}
#'   \item{config}{Screen method, phenotype class, label type, and parameters}
#'   \item{description}{Additional user-provided description}
#' }
#'
#' @details
#' Only for SigBridgeR version >= 3.7.0
#'
#'
#' The JSON file includes a header comment indicating it is required for cache identification and should not be modified.
#'
#' @examples
#' \dontrun{
#' WriteCacheMeta(
#'   file = "cache/meta.json",
#'   screen_method = "DEGAS",
#'   phenotype_class = "binary",
#'   params = list(n_top = 100),
#'   additional_description = "Test run"
#' )
#' }
#' @family cache_config
#' @name WriteCacheMeta
#' @export
WriteCacheMeta <- new_generic(
  name = "WriteCacheMeta",
  dispatch_args = "cache_config"
)

#' @rdname WriteCacheMeta
#' @export
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

#' @rdname WriteCacheMeta
#' @export
method(
  WriteCacheMeta,
  ScreenMethodConfig
) <- WriteCacheMeta.ScreenMethodConfig <- function(
  cache_config,
  file,
  additional_description = NULL,
  verbose = TRUE,
  ...
) {
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

#' @rdname WriteCacheMeta
#' @export
method(WriteCacheMeta, ScreenMethodCache) <- function(
  cache_config,
  additional_description = NULL,
  verbose = TRUE,
  ...
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
