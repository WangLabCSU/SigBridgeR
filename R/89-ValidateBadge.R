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
    if (value == 0L) {
      col <- cli::col_green
      sym <- cli_symbol$tick
      label <- name
    } else {
      col <- if (name == "note") cli::col_blue else cli::col_red
      sym <- cli_symbol$cross
      label <- if (value == 1L) name else paste0(name, "s")
    }

    col(paste(value, label, sym))
  }

  symbol_error <- symbol(check_res$error, "error")
  symbol_warn <- symbol(check_res$warn, "warning")
  symbol_note <- symbol(check_res$note, "note")

  cli::cli_text(symbol_error, " | ", symbol_warn, " | ", symbol_note)
}
