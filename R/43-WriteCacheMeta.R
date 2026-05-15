#' @title Write Cache Metadata
#'
#' @description Writes metadata for SigBridgeR cache identification to a JSON file.
#'
#' @param file Character string. Path to the JSON metadata file to be written.
#' @param screen_method Character string. Screening method name. Must match one of the keys in `ScreenStrategy`.
#' @param phenotype_class Character string. Type of phenotype data. One of `"binary"`, `"survival"`, or `"continuous"`.
#' @param label_type Character string. Type of label for the screening. Defaults to `screen_method`.
#' @param params List. Parameters used for the screening configuration.
#' @param additional_description Character string. Optional additional description text to include in the metadata. Default is `NULL`.
#' @param verbose Logical. Whether to print metadata preview for debugging. Default is `FALSE`.
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
#' @export
WriteCacheMeta <- function(
  file = "cache_config.json",
  screen_method,
  phenotype_class = c("binary", "survival", "continuous"),
  label_type = screen_method,
  params,
  additional_description = NULL,
  verbose = FALSE,
  ...
) {
  rlang::check_dots_empty0()
  rlang::check_installed("jsonlite", "writting .json config")
  chk::chk_string(file)
  chk::chk_string(phenotype_class)
  chk::chk_string(label_type)
  chk::chk_list(params)
  if (!is.null(additional_description)) {
    chk::chk_string(additional_description)
  }
  chk::chk_flag(verbose)
  chk::chk_length(file)

  screen_method <- SigBridgeRUtils::MatchArg(
    screen_method,
    names(ScreenStrategy),
    NULL
  )
  method_config <- ScreenStrategy[[screen_method]]

  phenotype_class <- SigBridgeRUtils::MatchArg(
    phenotype_class,
    c("binary", "survival", "continuous"),
    NULL
  )
  chk::chk_length(screen_method)
  chk::chk_length(label_type)

  cache_meta <- list(
    general = rlang::list2(
      file = file.path(getwd(), file),
      os = get_os_info(),
      time = get_time(),
      !!!get_r_info(),
      SigBridgeR_version = get_pkg_version()
    ),
    config = list(
      screen_method = screen_method,
      phenotype_class = phenotype_class,
      label_type = label_type,
      params = params
    ),
    description = additional_description
  )

  if (verbose) {
    cli::cli_h2(cli::col_blue("SigBridgeR Cache Metadata Preview"))
    cli::cli_dl(cache_meta)
    cli::cli_alert_info("Meta to: {.path {file}}")
  }

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

  json_content <- jsonlite::toJSON(cache_meta, auto_unbox = TRUE, pretty = TRUE)

  con <- file(file, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(json_content, con)

  if (verbose) {
    cli::cli_alert_success("Successfully wrote cache meta data.")
  }

  invisible(cache_meta)
}

#' @keywords internal
get_os_info <- function() {
  sys_info <- Sys.info()
  os_info <- trimws(paste(sys_info["sysname"], sys_info["release"]))
  if (!nzchar(os_info)) {
    return(sys_info["sysname"])
  }
  os_info
}

#' @keywords internal
get_r_info <- function() {
  list(
    r_version = R.version.string,
    arch = R.version$arch,
    platform = R.version$platform
  )
}

#' @keywords internal
get_time <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

#' @keywords internal
get_pkg_version <- function() {
  tryCatch(
    as.character(packageVersion("SigBridgeR")),
    error = function(e) {
      # devtools::load_all() 会将工作目录切换至包根目录
      read.dcf("DESCRIPTION")[, "Version"]
    }
  )
}

#' @keywords internal
get_cache_file_name <- function(
  screen_method,
  phenotype_class
) {
  paste(
    screen_method,
    phenotype_class,
    format(Sys.time(), "%Y%m%d%H%M"),
    sep = "_"
  )
}
