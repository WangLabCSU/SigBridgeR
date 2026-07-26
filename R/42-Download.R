#' @keywords internal
#' @inheritParams curl::curl_download
fileDownload <- function(
  url,
  destfile,
  quiet = !SigBridgeRUtils::getFuncOption("verbose"),
  mode = "wb"
) {
  check_installed("processx")
  download_res <- processx::run(
    "wget",
    c(
      "-c",
      url,
      "-O",
      destfile
    ),
    wd = getwd(),
    timeout = SigBridgeRUtils::getFuncOption("timeout"),
    echo_cmd = quiet,
    echo = quiet,
  )
  if (download_res$status == 0) {
    return(invisible())
  }
  if (!quiet) {
    cli::cli_alert_danger(
      "Download failed, trying {.fun download.file}"
    )
  }
  res <- utils::download.file(
    url = url,
    destfile = destfile,
    mode = mode,
    quiet = quiet
  )

  if (res == 0) {
    return(invisible())
  }

  if (!quiet) {
    cli::cli_alert_danger(
      "Download failed, trying {.fun curl::curl_download}"
    )
  }

  check_installed("curl")
  curl::curl_download(
    url = url,
    destfile = destfile,
    quiet = quiet,
    mode = mode
    # ,handle = handle
  )
}
