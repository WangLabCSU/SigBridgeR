#' @keywords internal
validate_warn <- function(message, ...) {
  cli::cli_inform(
    cli::col_magenta(paste(
      cli::symbol$pointer,
      message,
      "... WARNING"
    )),
    ...
  )
}

#' @keywords internal
validate_error <- function(message, ...) {
  cli::cli_inform(
    cli::col_red(paste(
      cli::symbol$pointer,
      message,
      "... ERROR"
    )),
    ...
  )
}

#' @keywords internal
validate_note <- function(message, ...) {
  cli::cli_inform(
    cli::col_white(paste(
      cli::symbol$pointer,
      message,
      "... NOTE"
    )),
    ...
  )
}

#' @keywords internal
validate_explain <- function(message, ...) {
  for (msg in message) {
    cli::cli_inform(c(" " = msg), ...)
  }
}

#' @keywords internal
validate_explain_speaker <- function(
  message = vector(),
  message_name = character(),
  type = c("error", "warn", "note"),
  ...
) {
  type <- match.arg(type)
  color_fun <- switch(
    type,
    error = cli::col_red,
    warn = cli::col_magenta,
    note = cli::col_white
  )

  prefix <- paste0(color_fun(message_name), ": ")
  prefix_width <- cli::ansi_nchar(prefix)

  # 2. 处理单行消息
  if (length(message) <= 1) {
    cli::cli_text(paste0(prefix, message))
    return(invisible(NULL))
  }
  cli::cli_text(paste0(prefix, message[1]))

  cli::cli_div(
    id = "validate_explain_speaker",
    theme = list(div = list("margin-left" = prefix_width))
  )
  on.exit(cli::cli_end(id = "validate_explain_speaker"), add = TRUE)

  for (msg in message[-1]) {
    cli::cli_ul(msg)
  }
}

#' @keywords internal
validate_success <- function(text, ...) {
  cli::cli_alert_success(cli::col_grey(text), ...)
}
