#' @title Generate a Template Screening Function with Optional Roxygen2 Documentation
#' @family Add_Screen_method
#' @description
#' Creates a scaffold for a new screening function that conforms to the expected
#' interface for bulk-scRNA integration workflows. The template includes standard
#' parameters (e.g., `matched_bulk`, `sc_data`, `phenotype`) and a return structure
#' compatible with downstream validation (e.g., `ValidateScreenFunc`). Optionally
#' includes roxygen2 documentation and supports writing to a new or existing `.R` file.
#'
#' When \code{filename = "current"}, the function attempts to append to the active
#' script in RStudio or Positron (requires \code{rstudioapi}).
#'
#' @param filename Character. Path to the output `.R` file. If \code{"current"} (case-insensitive),
#'   uses the currently open script in RStudio/Positron. Must end in `.R` if specified explicitly.
#'   Default: \code{"my_screen.R"}.
#' @param func_name Character. Name of the generated function. Must be a valid R identifier
#'   (starts with letter, followed by letters, digits, `_`, or `.`). Default: \code{"my_screen"}.
#' @param append Logical. If \code{TRUE} (default), appends to an existing file; otherwise, overwrites it.
#'   If overwriting a non-empty file, prompts for user confirmation.
#' @param documentation Logical. Whether to include a full roxygen2 documentation block above the function.
#'   Default: \code{FALSE}. If \code{TRUE}, users should run \code{devtools::document()} afterward.
#' @param open Logical. Whether to open the file in the RStudio/Positron editor after writing (if available).
#'   Default: \code{TRUE}.
#' @param ... Additional arguments (currently unused, reserved for future extension).
#'
#' @return Invisibly returns the path to the written file (\code{filename}).
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a new screen function with documentation
#' TemplateScreenFunc(
#'   filename = "custom_screens.R",
#'   func_name = "MyCoxScreen",
#'   documentation = TRUE
#' )
#'
#' # Append to current script in RStudio
#' TemplateScreenFunc(filename = "current", func_name = "QuickFilter")
#' }
TemplateScreenFunc <- function(
  filename = NULL,
  func_name = "my_screen",
  append = TRUE,
  documentation = FALSE,
  open = is_interactive(),
  ...
) {
  filename <- filename %||% "my_screen.R"
  chk::chk_character(filename)
  chk::chk_length(filename) # 1
  chk::chk_not_null(filename)

  if (tolower(trimws(filename)) == "current") {
    # * Check rstudioapi availability
    rlang::check_installed("rstudioapi")

    if (!rstudioapi::isAvailable()) {
      Abort(
        "{.pkg rstudioapi} is not available. This function requires RStudio or Positron session when specified filename as `current`.",
        "Please specify a filename directly instead of `current`"
      )
    }

    # Get current script path
    context <- rstudioapi::getSourceEditorContext()
    filename <- context$path
    cli::cli_alert_info("Using current script: {.file {basename(filename)}}")
  } else {
    # * Validate filename is character and has .R extension
    if (!grepl("\\.R$", filename, ignore.case = TRUE)) {
      Abort(
        "`filename` must have a {.path .R} extension (case-insensitive)."
      )
    }
  }

  # * validate function name
  if (!grepl("^[a-zA-Z][a-zA-Z0-9._]*$", func_name)) {
    Abort(
      "Invalid `func_name`: {.val {func_name}}"
    )
  }

  # * Check file existence and content
  file_exists <- file.exists(filename)
  file_not_empty <- file_exists && (file.info(filename)$size > 0)

  # * Interactive confirmation for overwriting non-empty files
  if (!append && file_not_empty) {
    response <- utils::askYesNo("Overwrite existing content? ")

    if (!isTRUE(response)) {
      Abort("Operation cancelled by user.")
    }
  }

  # * Build template function body from template file
  template_body <- build_template_body(func_name, documentation)

  write_mode <- if (append && file_exists) "a" else "w"
  con <- file(filename, open = write_mode)
  on.exit(close(con), add = TRUE)

  # Add newline before appending to non-empty files for separation
  if (append && file_not_empty) {
    writeLines("\n", con = con)
  }

  writeLines(template_body, con = con, sep = "\n")

  # Success message
  action <- if (append && file_exists) "appended to" else "written to"
  cli::cli_alert_success(
    "Template function {.fun {func_name}} successfully {action} {.file {basename(filename)}}"
  )

  if (documentation) {
    cli::cli_alert_success(
      "Roxygen2 documentation included (use {.code devtools::document()} to generate NAMESPACE)."
    )
  }

  # Open file in RStudio if available
  if (open && rstudioapi::isAvailable()) {
    rstudioapi::navigateToFile(filename)
    cli::cli_text(
      "{cli::col_red(cli::symbol$checkbox_off)} File opened in editor."
    )
  }

  invisible(filename)
}

# Reads the template file and constructs the output with optional documentation.
# @param func_name Character. The function name to substitute into the template.
# @param documentation Logical. Whether to include roxygen2 documentation.
# @return Character vector of lines.
# @keywords internal
build_template_body <- function(func_name, documentation) {
  template_path <- system.file(
    "template",
    "template_screen_fun.R",
    package = "SigBridgeR",
    mustWork = TRUE
  )

  lines <- readLines(template_path, warn = FALSE)

  # Substitute function name placeholder
  lines <- gsub("FUNC_NAME", func_name, lines, fixed = TRUE)

  if (!documentation) {
    # Remove roxygen2 documentation lines (lines starting with #')
    lines <- lines[!grepl("^#'", lines)]
  }

  lines
}
