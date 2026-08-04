#' @keywords internal
#' @inheritParams curl::curl_download
fileDownload <- function(
  url,
  destfile,
  quiet = !SigBridgeRUtils::getFuncOption("verbose"),
  mode = "wb",
  retries = 3L,
  timeout = 600
) {
  check_installed("curl")
  dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)

  tmp <- paste0(destfile, ".tmp")
  on.exit(
    {
      if (file.exists(tmp)) unlink(tmp)
    },
    add = TRUE
  )

  h <- curl::new_handle(
    followlocation = TRUE,
    connecttimeout = 20,
    timeout = timeout,
    useragent = sprintf(
      "R/%s R (%s)",
      getRversion(),
      paste(
        getRversion(),
        R.version["platform"],
        R.version["arch"],
        R.version["os"]
      )
    )
  )

  last_err <- NULL

  for (i in seq_len(retries)) {
    if (file.exists(tmp)) {
      unlink(tmp)
    }

    err <- tryCatch(
      {
        curl::curl_download(
          url = url,
          destfile = tmp,
          quiet = quiet,
          mode = mode,
          handle = h
        )
        NULL
      },
      error = function(e) e
    )

    if (is.null(err)) {
      ok <- file.rename(tmp, destfile)

      if (!ok) {
        ok <- file.copy(tmp, destfile, overwrite = TRUE)
        unlink(tmp)
      }

      if (ok) {
        return(invisible(destfile))
      }

      last_err <- simpleError(
        sprintf("Failed to move temporary file to %s", destfile)
      )
    } else {
      last_err <- err
    }

    if (!quiet) {
      cli::cli_alert_warning(
        "Download attempt {i}/{retries} failed: {conditionMessage(last_err)}"
      )
    }

    Sys.sleep(min(2^i, 10))
  }

  stop(
    "Download failed after ",
    retries,
    " attempts: ",
    conditionMessage(last_err),
    call. = FALSE
  )
}
