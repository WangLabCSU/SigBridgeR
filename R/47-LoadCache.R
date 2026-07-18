#' @title Load a Cached R Object
#'
#' @description
#' Reads a cached R object from disk. The file format is determined by the
#' file extension:
#' \itemize{
#'   \item \code{.qs2}: loaded via \code{qs2::qs_read()}
#'   \item \code{.RData} / \code{.rds}: loaded via \code{readRDS()}
#'   \item \code{.csv}: loaded via \code{data.table::fread()}
#' }
#'
#' If the file was originally saved in \code{.qs2} format but \code{qs2} is
#' not installed on the current machine, a clear error message is shown with
#' installation instructions.
#'
#' @param file Character string. Path to the cache file to load.
#' @param ... Additional arguments (must be empty, checked by
#'   \code{rlang::check_dots_empty0()}).
#'
#' @return The cached R object.
#'
#' @family cache_config
#' @export
#'
#' @examples
#' \dontrun{
#' result <- LoadCache("Scissor_res/survival_2025_01_01/result.qs2")
#' feature_df <- LoadCache("Scissor_res/survival_2025_01_01/feature_table.csv")
#' }
LoadCache <- function(file, ...) {
  rlang::check_dots_empty0()
  chk::chk_string(file)
  chk::chk_file(file)

  ext <- tools::file_ext(file)

  result <- switch(
    ext,
    "qs2" = {
      if (!rlang::is_installed("qs2")) {
        cli::cli_abort(c(
          "x" = "Package {.pkg qs2} is required to read this cache file.",
          ">" = "Install it with {.code install.packages('qs2')}.",
          "i" = "This file was previously cached in qs2 format."
        ))
      }
      qs2::qs_read(file)
    },
    RData = readRDS(file),
    rds = readRDS(file),
    csv = data.table::fread(file),
    cli::cli_abort(c(
      "x" = "Unsupported file extension: {.val {ext}}",
      ">" = "Supported extensions: {.val {c('qs2', 'RData', 'rds', 'csv')}}"
    ))
  )

  result
}
