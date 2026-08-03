#' @title Write an R Object to Cache
#'
#' @description
#' Writes an R object to the cache directory. Uses `qs2` for fast
#' serialization when available, falling back to `rdata` for general
#' objects. Small data frames and data.tables are saved as CSV via
#' `data.table::fwrite()` for easy inspection.
#'
#' # Methods
#' `WriteCache` is an S7 generic with methods available for the following
#' classes:
#'
#' `r doclisting::methods_list("WriteCache")`
#'
#' @details
#' The `ScreenMethodCache` method writes both the cache metadata (via
#' [WriteCacheMeta()]) and the cached data objects to disk. The
#' `ScreenMethodConfig` method only writes the metadata.
#'
#' ## Format detection (ScreenMethodCache)
#'
#' When `format = "auto"` (default), the format for each data object is
#' detected automatically:
#'
#' * If the object is a 2D structure (data.frame/matrix) with at most
#'   `max_rows` rows and `max_cols` columns, CSV format is used.
#' * If `qs2` is installed and the object is a list, vector,
#'   data.frame, or matrix, `qdata` format is used.
#' * If `qs2` is installed but the object is not qdata-compatible,
#'   `qs2` format is used.
#' * Otherwise, `rdata` format is used as a fallback.
#'
#' @param cache A `ScreenMethodConfig` or `ScreenMethodCache` S7 object
#'   containing the data and configuration to write.
#' @param ... Additional arguments passed to methods.
#'
#'   For `ScreenMethodConfig`:
#'   - `file`: `character(1)`. Path to the output `cache_config.json` file.
#'   - `additional_description`: `character` or `NULL`. An optional string
#'     describing the cache entry. Default: `NULL`.
#'   - `verbose`: `logical`. Whether to print progress messages. Default:
#'     `TRUE`.
#'
#'   For `ScreenMethodCache`:
#'   - `format`: `character`. Cache format. One of `"auto"`, `"qs2"`,
#'     `"qdata"`, `"rdata"`, or `"csv"`. When `"auto"` default, the format is
#'     detected per object.
#'   - `max_rows`: `integer`. Maximum number of rows to consider a data.frame
#'     small for CSV output. Default: `1000L`.
#'   - `max_cols`: `integer`. Maximum number of columns to consider a
#'     data.frame small for CSV output. Default: `20L`.
#'   - `additional_description`: `character` or `NULL`. An optional string
#'     describing the cache entry. Default: `NULL`.
#'   - `verbose`: `logical`. Whether to print progress messages. Default:
#'     `TRUE`.
#'
#' @returns Invisible. Returns the absolute path(s) to the written cache
#'   file(s).
#'
#' @family cache_config
#' @name WriteCache
#' @export
#'
#' @examples
#' \dontrun{
#' # Write using a ScreenMethodCache
#' cache <- ScreenMethodCache(
#'   cache_path = "Scissor_res/survival_2025_01_01/",
#'   cache_config_path = "Scissor_res/survival_2025_01_01/cache_config.json",
#'   cache_data = list(result = my_matrix),
#'   screen_method_config = config
#' )
#' WriteCache(cache)
#'
#' # Force CSV format for small data frames
#' WriteCache(cache, format = "csv")
#' }
WriteCache <- new_generic(
  "WriteCache",
  dispatch_args = "cache"
)

method(WriteCache, class_any) <- function(
  cache,
  ...
) {
  cls_cache <- class(cache)
  expected_cls <- c("ScreenMethodConfig", "ScreenMethodCache")
  Abort(
    "cache must be a class of {.cls {expected_cls}}",
    "Current class is {.cls {cls_cache}}"
  )
}

method(WriteCache, ScreenMethodConfig) <- function(
  cache,
  ...,
  additional_description = NULL,
  verbose = TRUE
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
  ...,
  format = c("auto", "qs2", "qdata", "rdata", "csv"),
  max_rows = 1000L,
  max_cols = 20L,
  additional_description = NULL,
  verbose = TRUE
) {
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
  saving_fn <- purrr::in_parallel(function(x, file, format) {
    WriteSingleCacheData(
      x = x,
      file = file,
      format = format
    )
  })
  purrr::pwalk(
    list(cache_data, files, formats),
    saving_fn,
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
