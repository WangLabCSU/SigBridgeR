# nocov start

EnsureParentDir <- function(path) {
  parent <- dirname(normalizePath(path, mustWork = FALSE))

  if (dir.exists(parent)) {
    return(invisible(TRUE))
  }

  ans <- utils::askYesNo(
    cli::cli_fmt(cli::cli_inform(
      "{cli::symbol$menu} Create parent directory {.file {parent}}?"
    ))
  )

  if (!tolower(ans) %in% c("y", "yes")) {
    Abort("Canceled by user")
  }

  dir.create(parent, recursive = TRUE, showWarnings = FALSE)

  invisible(TRUE)
}
