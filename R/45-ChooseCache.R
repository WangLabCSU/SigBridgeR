#' @title Choose a Cache Directory
#'
#' @description
#' Selects a cache directory from the specified path. If only one cache
#' subdirectory exists, it is returned automatically; if multiple are found,
#' an interactive menu is presented for the user to choose.
#'
#' @details
#' The function expects the parent directory name to follow the convention
#' `{method_name}_res` (e.g., `"Scissor_res"`, `"scPAS_res"`). It validates
#' the directory name against the registered screening methods in
#' `ScreenStrategy`.
#'
#' If no cache subdirectories are found, the function aborts with an error.
#' If exactly one subdirectory exists, it is returned without prompting.
#' If multiple subdirectories exist, `utils::menu()` is used to present an
#' interactive selection menu.
#'
#' @param directory `character`. Path to the parent directory containing
#'   cache subdirectories. The directory name must match
#'   `{method_name}_res` for a registered screening method.
#'
#' @returns A `character` string: the path to the selected cache directory.
#'
#' @family cache_config
#' @export
#'
#' @examples
#' \dontrun{
#' cache_path <- ChooseCache("./Scissor_res")
#' }
ChooseCache <- function(directory) {
  expected_dir_name <- paste0(names(ScreenStrategy), "_res")
  if (!basename(directory) %chin% expected_dir_name) {
    Abort(
      "{.path {directory}} is not a cache directory",
      "Expected directory name is {.file {expected_dir_name}}"
    )
  }

  if (!dir.exists(directory)) {
    Abort("{.path {directory}} not exists")
  }
  cache_dirs <- list.dirs(directory, recursive = FALSE)
  if (length(cache_dirs) == 0) {
    Abort("no cache found in {.path {directory}}")
  } else if (length(cache_dirs) == 1L) {
    return(cache_dirs[[1L]])
  }

  # length > 1
  choice <- utils::menu(
    cache_dirs,
    title = cli::cli_fmt(cli::cli_inform(
      "{cli::symbol$menu} Which cache do you want to use?"
    ))
  )

  n_cache_dirs <- length(cache_dirs)
  if (!is.numeric(choice) || choice < 1 || choice > n_cache_dirs) {
    Abort(
      "Invalid choice: {.val {choice}}",
      "Expected a number between 1 and {n_cache_dirs}"
    )
  }

  cache_dirs[[choice]]
}
