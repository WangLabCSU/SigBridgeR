ValidateBadge <- function(...) {
  merge_check_res <- function(...) {
    lists <- list2(...)

    all_names <- c("error", "warn", "note")
    lists <- lapply(lists, function(x) {
      missing <- setdiff(all_names, names(x))
      if (length(missing)) {
        x[missing] <- 0L
      }
      x[all_names]
    })

    purrr::reduce(lists, ~ purrr::map2(.x, .y, `+`))
  }
  check_res <- merge_check_res(...)

  symbol <- function(value = integer(), name = character()) {
    t <- cli::symbol$tick
    x <- cli::symbol$cross

    if (value == 0L) {
      cli::cli_fmt(cli::cli_text(cli::col_green("{value} {name} {t}")))
    } else if (name == "note") {
      if (value == 1L) {
        cli::cli_fmt(cli::cli_text(cli::col_blue("{value} {name} {x}")))
      } else {
        cli::cli_fmt(cli::cli_text(cli::col_blue("{value} {name}s {x}")))
      }
    } else if (value == 1L) {
      cli::cli_fmt(cli::cli_text(cli::col_red("{value} {name} {x}")))
    } else {
      cli::cli_fmt(cli::cli_text(cli::col_red("{value} {name}s {x}")))
    }
  }

  symbol_error <- symbol(check_res$error, "error")
  symbol_warn <- symbol(check_res$warn, "warning")
  symbol_note <- symbol(check_res$note, "note")

  cli::cli_div()
  cli::cli_text(symbol_error, " | ", symbol_warn, " | ", symbol_note)
  cli::cli_end()
}
