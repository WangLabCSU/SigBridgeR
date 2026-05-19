#' @title Write an R Object to Cache
#'
#' @description
#' Writes an R object to the cache directory. Uses \code{qs2} for fast
#' serialization when available, falling back to \code{RDS} for general
#' objects. Small data frames and data.tables are saved as CSV via
#' \code{data.table::fwrite()} for easy inspection.
#'
#' @param x An R object to be cached.
#' @param file Character string. Path to the cache file (with or without
#'   extension). If no extension is given, it is appended automatically.
#' @param format Character string. Cache format. One of \code{"auto"},
#'   \code{"qs2"}, \code{"rds"}, or \code{"csv"}. When \code{"auto"}
#'   (default):
#'   \itemize{
#'     \item If \code{x} is a data.frame with at most \code{max_rows} rows
#'       and \code{max_cols} columns, CSV format is used.
#'     \item Otherwise, \code{qs2} is preferred when the package is
#'       installed, falling back to \code{RDS}.
#'   }
#' @param max_rows Integer. Maximum number of rows to consider a data.frame
#'   "small". Default is \code{1000L}.
#' @param max_cols Integer. Maximum number of columns to consider a
#'   data.frame "small". Default is \code{20L}.
#' @param ... Additional arguments (must be empty, checked by
#'   \code{rlang::check_dots_empty0()}).
#'
#' @return Invisible. Returns the absolute path to the written cache file.
#'
#' @family cache_config
#' @export
#'
#' @examples
#' \dontrun{
#' # Write a matrix as qs2
#' WriteCache(cache_res, "Scissor_res/survival_2025_01_01/result.qs2")
#'
#' # Write a small data.frame as CSV
#' WriteCache(small_df, "Scissor_res/survival_2025_01_01/feature_table.csv")
#' }
WriteCache <- function(
  x,
  file,
  format = c("auto", "qs2", "rds", "csv"),
  max_rows = 1000L,
  max_cols = 20L,
  ...
) {
  rlang::check_dots_empty0()
  chk::chk_string(file)
  chk::chk_count(max_rows)
  chk::chk_count(max_cols)

  format <- SigBridgeRUtils::MatchArg(format, c("auto", "qs2", "rds", "csv"))

  # -- auto-detect format ---------------------------------------------------
  if (format == "auto") {
    is_small_df <- is.data.frame(x) &&
      NROW(x) <= max_rows &&
      NCOL(x) <= max_cols
    if (is_small_df) {
      format <- "csv"
    } else if (rlang::is_installed("qs2")) {
      format <- "qs2"
    } else {
      format <- "rds"
    }
  }

  # -- normalise extension --------------------------------------------------
  ext <- tools::file_ext(file)
  if (!nzchar(ext)) {
    file <- switch(
      format,
      "qs2" = paste0(file, ".qs2"),
      "rds" = paste0(file, ".RData"),
      "csv" = paste0(file, ".csv")
    )
  }

  # -- ensure parent directory exists ---------------------------------------
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

  # -- write ----------------------------------------------------------------
  switch(
    format,
    "qs2" = {
      rlang::check_installed(
        pkg = "qs2",
        reason = "to write cache in qs2 format"
      )
      qs2::qs_save(x, file)
    },
    "rds" = {
      saveRDS(object = x, file = file)
    },
    "csv" = {
      data.table::fwrite(x = x, file = file)
    }
  )

  invisible(normalizePath(file, mustWork = FALSE))
}
