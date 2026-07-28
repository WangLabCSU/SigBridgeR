#' @title Write an R Object to Cache
#'
#' @description
#' Writes an R object to the cache directory. Uses \code{qs2} for fast
#' serialization when available, falling back to \code{rdata} for general
#' objects. Small data frames and data.tables are saved as CSV via
#' \code{data.table::fwrite()} for easy inspection.
#'
#' @param x An R object to be cached.
#' @param file Character string. Path to the cache file (with or without
#'   extension). If no extension is given, it is appended automatically.
#' @param format Character string. Cache format. One of \code{"auto"},
#'   \code{"qs2"}, \code{"rdata"}, or \code{"csv"}. When \code{"auto"}
#'   (default):
#'   \itemize{
#'     \item If \code{x} is a data.frame with at most \code{max_rows} rows
#'       and \code{max_cols} columns, CSV format is used.
#'     \item Otherwise, \code{qs2} is preferred when the package is
#'       installed, falling back to \code{rdata}.
#'   }
#' @param max_rows Integer. Maximum number of rows to consider a data.frame
#'   "small". Default is \code{1000L}.
#' @param max_cols Integer. Maximum number of columns to consider a
#'   data.frame "small". Default is \code{20L}.
#' @param ... Additional arguments (must be empty, checked by
#'   \code{check_dots_empty0()}).
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
WriteCache <- new_generic("WriteCache", dispatch_args = "cache")

method(WriteCache, ScreenMethodConfig) <- function(
  cache,
  additional_description = NULL,
  verbose = TRUE,
  ...
) {
  WriteCacheMeta(
    cache_config = cache@screen_method_config,
    file = cache@cache_config_path,
    additional_description = additional_description,
    verbose = verbose,
    ...
  )
}

method(WriteCache, ScreenMethodCache) <- function(
  cache,
  format = c("auto", "qs2", "qdata", "rdata", "csv"),
  max_rows = 1000L,
  max_cols = 20L,
  additional_description = NULL,
  verbose = TRUE,
  ...
) {
  check_dots_empty()

  cache_path <- cache@cache_path
  cache_data <- cache@cache_data

  chk::chk_chr(cache_path)
  chk::chk_count(max_rows)
  chk::chk_count(max_cols)

  # -- write cache meta -------------------------------------------
  WriteCacheMeta(
    cache_config = cache@screen_method_config,
    file = cache@cache_config_path,
    additional_description = additional_description,
    verbose = verbose,
    ...
  )

  format <- SigBridgeRUtils::MatchArg(
    format,
    c("auto", "qs2", "qdata", "rdata", "csv")
  )

  # -- detect format --------------------------------------------------------

  formats <- if (identical(format, "auto")) {
    purrr::map_chr(cache_data, \(x) detect_format(x, max_rows, max_cols))
  } else {
    rep(format, length(cache_data))
  }

  if (any(c("qs2", "qdata") %chin% formats)) {
    check_installed(
      pkg = "qs2",
      reason = "to write cache in qs2/qdata format"
    )
  }

  # -- cache_path as directory ----------------------------------------------
  cache_dir <- if (nzchar(tools::file_ext(cache_path))) {
    tools::file_path_sans_ext(cache_path)
  } else {
    cache_path
  }

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  # -- file names -----------------------------------------------------------
  cache_names <- names(cache_data)

  files <- file.path(
    cache_dir,
    paste0(cache_names, formats)
  )

  # -- batch write ----------------------------------------------------------
  purrr::pwalk(
    list(cache_data, files, formats),
    purrr::in_parallel(function(x, file, format) {
      WriteSingleCacheData(
        x = x,
        file = file,
        format = format
      )
    }),
    .progress = if (isTRUE(verbose)) "Saving cache files" else FALSE
  )

  invisible(files)
}

detect_format <- function(x, max_rows, max_cols) {
  is_2d_obj <- !is.null(dim(x)) && length(dim(x)) == 2L

  is_small_df <- is_2d_obj &&
    NROW(x) <= max_rows &&
    NCOL(x) <= max_cols

  is_qdata_allowed <- is.list(x) ||
    is.vector(x) ||
    is.data.frame(x) ||
    is.matrix(x)

  if (is_small_df) {
    "csv"
  } else if (is_installed("qs2") && is_qdata_allowed) {
    "qdata"
  } else if (is_installed("qs2")) {
    "qs2"
  } else {
    "rdata"
  }
}

WriteSingleCacheData <- function(x, file, format) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

  switch(
    format,

    "qs2" = ,
    "qdata" = {
      qs2::qs_save(x, file)
    },

    "rdata" = {
      saveRDS(object = x, file = file)
    },

    "csv" = {
      x <- data.table::as.data.table(x)
      data.table::fwrite(x = x, file = file)
    }
  )

  invisible(file)
}
