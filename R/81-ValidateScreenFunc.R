#' @title Validate Custom Screening Function Compliance
#' @family Add_Screen_method
#'
#' @description
#' Verifies if a user-provided function meets the interface requirements for
#' the screening pipeline.
#'
#' @param func Function to validate. Must accept core parameters and return
#'   expected structure.
#' @param ... For future updates
#'
#' @return \code{TRUE} if all checks pass, otherwise terminates with diagnostic report
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example compliant function
#' my_screen <- function(
#'   sc_data,
#'   matched_bulk,
#'   phenotype,
#'   label_type,
#'   phenotype_class,
#'   verbose = FALSE,
#'   ...
#' ) {
#'   if (verbose) {
#'     message("Running custom screen")
#'   }
#'   list(
#'     scRNA_data = sc_data, # Must be Seurat object
#'     results = data.frame(score = runif(10))
#'   )
#' }
#'
#' ValidateScreenFunc(my_screen)
#' }
ValidateScreenFunc <- function(
  func,
  ...
) {
  cli::cli_h1(cli::col_cyan("Screening Function Validation"))

  if (is_installed("tictoc")) {
    tictoc::tic()
  }

  cli::cli_text(sprintf(
    "Start at %s",
    format(Sys.time(), "%Y/%m/%d %H:%M:%S")
  ))
  cli::cat_line()

  inputs_check <- ValidateArgsInputs(func, ...)
  cli::cat_line()

  verbose_check <- ValidateArgsVerbose(func, ...)
  cli::cat_line()

  syntax_check <- ValidateArgsSyntax(func, ...)
  cli::cat_line()

  return_check <- ValidateReturn(func, ...)
  cli::cat_line()

  dir_name_check <- ValidateDirName(func, ...)
  cli::cat_line()

  if (is_installed("tictoc")) {
    end_time <- tictoc::toc(quiet = TRUE)
    cli::cli_text(tictoc::toc.outmsg(end_time$tic, end_time$toc, "Duration"))
    cli::cat_line()
  }

  ValidateBadge(
    inputs_check,
    verbose_check,
    syntax_check,
    return_check,
    dir_name_check
  )
  cli::cat_line()

  invisible(TRUE)
}
