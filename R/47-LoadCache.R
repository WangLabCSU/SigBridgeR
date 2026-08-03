#' @title Load a Cached R Object
#'
#' @description
#' Reads a cached R object from disk. The file format is determined by the
#' file extension:
#'
#' * `.qs2`: loaded via `qs2::qs_read()`
#' * `.RData` / `.rds`: loaded via `readRDS()`
#' * `.csv`: loaded via `data.table::fread()`
#'
#' If the file was originally saved in \code{.qs2} format but \code{qs2} is
#' not installed on the current machine, a clear error message is shown with
#' installation instructions.
#'
#' @param file Character string. Path to the cache file to load.
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
LoadCache <- function(file) {
  chk::chk_string(file)
  chk::chk_file(file)

  ext <- tools::file_ext(file)

  switch(
    ext,
    "qs2" = ,
    "qdata" = {
      check_installed("qs2")
      qs2::qs_read(file)
    },
    RData = readRDS(file),
    rds = readRDS(file),
    csv = data.table::fread(file),
    Abort(
      "Unsupported file extension: {.val {ext}}",
      "Supported extensions: {.val {c('qs2', 'qdata', 'RData', 'rds', 'csv')}}"
    )
  )
}
