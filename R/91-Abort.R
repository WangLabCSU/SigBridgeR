# nocov start

Abort <- function(
  error_msg,
  tips = NULL,
  info = NULL,
  type = NULL,
  .envir = caller_env()
) {
  if (!is.null(type)) {
    error_msg <- paste0(cli::col_red(type), " ", error_msg)
  }

  cli::cli_abort(
    message = c(
      "x" = error_msg,
      ">" = tips,
      "i" = info
    ),
    call = .envir,
    .envir = .envir,
    .frame = .envir
  )
}
