#' @title Download & Load Reference Data
#'
#' @description
#' This function checks if the data already exists in a local cache.
#' If not, it downloads the file from the remote repository using multiple sources with fallback.
#'
#' @param data_type The type of data to download. Must be one of "survival", "binary", or "continuous".
#' @param path Optional path to save the downloaded file, default: NULL, saving in package.
#' @param timeout Integer. Connection timeout in seconds (default: 60)
#' @param ... Additional arguments. Currently supports:
#'    - `verbose`: Logical indicating whether to print progress messages. Defaults to `TRUE`.
#'
#' @return The requested datasets, stored in a list.
#' @export
#'
#' @examples
#' \dontrun{
#' ref_data <- LoadRefData(data_type = c("survival"))
#' }
#'
LoadRefData <- function(
  data_type = c("survival", "binary", "continuous"),
  path = NULL,
  timeout = SigBridgeRUtils::getFuncOption("timeout"),
  ...
) {
  dots <- list2(...)
  verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption("verbose")
  chk::chk_whole_number(timeout)

  data_type <- SigBridgeRUtils::MatchArg(
    data_type,
    c("survival", "binary", "continuous"),
    NULL
  )
  if (!is.null(path)) {
    chk::chk_dir(path)
  } else {
    path <- if (.Platform$OS.type == "windows") {
      file.path(Sys.getenv("LOCALAPPDATA"), "SigBridgeR", "cache")
    } else {
      file.path(Sys.getenv("HOME"), ".cache", "SigBridgeR")
    }
    cli::cli_alert_info("Using default cache path: {.path {path}}")
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }

  local_file <- file.path(path, glue::glue("{data_type}_ref_data.rds"))

  if (!file.exists(local_file)) {
    if (verbose) {
      cli::cli_alert_info("Downloading reference data...")
    }

    # Define multiple sources with priority order
    data_urls <- list(
      survival = list(
        "GitHub Raw" = "https://raw.githubusercontent.com/WangLabCSU/SigBridgeR/refs/heads/main/inst/extdata/survival_example_data.rds",
        "GitHub API" = "https://raw.githubusercontent.com/WangLabCSU/SigBridgeR/main/inst/extdata/survival_example_data.rds",
        "CDN (jsDelivr)" = "https://cdn.jsdelivr.net/gh/WangLabCSU/SigBridgeR@main/inst/extdata/survival_example_data.rds"
      ),
      binary = list(
        "GitHub Raw" = "https://raw.githubusercontent.com/WangLabCSU/SigBridgeR/refs/heads/main/inst/extdata/binary_example_data.rds",
        "GitHub API" = "https://raw.githubusercontent.com/WangLabCSU/SigBridgeR/main/inst/extdata/binary_example_data.rds",
        "CDN (jsDelivr)" = "https://cdn.jsdelivr.net/gh/WangLabCSU/SigBridgeR@main/inst/extdata/binary_example_data.rds"
      ),
      continuous = list(
        "GitHub Raw" = "https://raw.githubusercontent.com/WangLabCSU/SigBridgeR/refs/heads/main/inst/extdata/continuous_example_data.rds",
        "GitHub API" = "https://raw.githubusercontent.com/WangLabCSU/SigBridgeR/main/inst/extdata/continuous_example_data.rds",
        "CDN (jsDelivr)" = "https://cdn.jsdelivr.net/gh/WangLabCSU/SigBridgeR@main/inst/extdata/continuous_example_data.rds"
      )
    )

    available_urls <- data_urls[[data_type]]
    success <- FALSE

    # Try each source in priority order
    for (i in seq_along(available_urls)) {
      source_name <- names(available_urls)[i]
      data_url <- available_urls[[i]]

      if (verbose) {
        cli::cli_alert_info("Trying {source_name}...")
      }

      success <- try_fetch(
        {
          fileDownload(
            url = data_url,
            destfile = local_file,
            mode = "wb",
            quiet = FALSE # Show progress
          )
          TRUE
        },
        error = function(e) {
          if (file.exists(local_file)) {
            unlink(local_file)
          }

          # If this was the last source, show final error
          if (i == length(available_urls)) {
            Abort(
              cli::col_red("All download attempts failed."),
              "Please check your internet connection or try again later.",
              "Error from last attempt: {e$message}"
            )
          }
          FALSE
        }
      )

      if (success) {
        if (verbose) {
          cli::cli_alert_success(cli::col_green("Data loaded successfully."))
        }
        break
      }
    }
  } else if (verbose) {
    cli::cli_alert_info("Found cached data ({.file {local_file}}).")
  }

  readRDS(local_file)
}
