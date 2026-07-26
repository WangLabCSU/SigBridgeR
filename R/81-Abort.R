Abort <- function(error_msg, tips = NULL, info = NULL, .envir = caller_env()) {
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
