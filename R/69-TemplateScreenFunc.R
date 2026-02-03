#' @title Generate a Template Screening Function with Optional Roxygen2 Documentation
#' @keywords Add_Screen_method
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
  filename = "my_screen.R",
  func_name = "my_screen",
  append = TRUE,
  documentation = FALSE,
  open = TRUE,
  ...
) {
  chk::chk_character(filename)
  chk::chk_length(filename) # 1

  if (is.null(filename)) {
    filename <- "my_screen.R"
    cli::cli_alert_info("Using default filename: {.file {filename}}")
  } else if (tolower(trimws(filename)) == "current") {
    # * Check rstudioapi availability
    rlang::check_installed("rstudioapi")

    if (!rstudioapi::isAvailable()) {
      cli::cli_abort(c(
        "x" = "{.pkg rstudioapi} is not available. This function requires RStudio or Positron session when specified filename as `current`.",
        ">" = "Please specify a filename directly instead of `current`"
      ))
    }

    # Get current script path
    context <- rstudioapi::getSourceEditorContext()
    filename <- context$path
    cli::cli_alert_info("Using current script: {.file {basename(filename)}}")
  } else {
    # * Validate filename is character and has .R extension
    if (!grepl("\\.R$", filename, ignore.case = TRUE)) {
      cli::cli_abort(
        "`filename` must have a {.path .R} extension (case-insensitive)."
      )
    }
  }

  # * validate function name
  if (!grepl("^[a-zA-Z][a-zA-Z0-9._]*$", func_name)) {
    cli::cli_abort(
      "Invalid `func_name`: {.val {func_name}}"
    )
  }

  # * Check file existence and content
  file_exists <- file.exists(filename)
  file_not_empty <- file_exists && (file.info(filename)$size > 0)

  # * Interactive confirmation for overwriting non-empty files
  if (!append && file_not_empty) {
    cli::cli_inform("Overwrite existing content? [y/N]: ", .width = "auto")

    response <- readline(prompt = "")
    if (tolower(response) != "y" && tolower(response) != "yes") {
      cli::cli_inform("Operation cancelled by user.")
      return(NULL)
    }
  }

  # * Build roxygen2 documentation if requested
  func_doc <- roxygen2_doc(documentation = documentation)

  # * Build template function body
  template_body <- func_body(func_doc = func_doc, func_name = func_name)

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

#' @keywords internal
roxygen2_doc <- function(documentation = FALSE) {
  if (!documentation) {
    return("")
  }

  paste(
    "#' @title Screen Function Template",
    "#'",
    "#' @param matched_bulk Matrix or data frame of preprocessed bulk RNA-seq expression",
    "#'        data (genes x samples). Column names must match names/IDs in `phenotype`.",
    "#' @param sc_data A matrix/Matrix (genes x cells) or a Seurat object containing scRNA-seq data to be screened.",
    "#' @param phenotype Phenotype data, either:",
    "#'        - Named vector (names match `matched_bulk` columns), or",
    "#'        - Patient survival Data frame with row names matching `matched_bulk` columns, colnames named \"time\" and \"status\"",
    "#' @param label_type Character specifying phenotype label type",
    "#' @param phenotype_class Type of phenotypic outcome (must be consistent with input data):",
    "#'        - `\"binary\"`: Binary traits (e.g., case/control)",
    "#'        - `\"continuous\"`: Continuous measurements",
    "#'        - `\"survival\"`: Survival infomation",
    "#' @param ... Additional arguments passed to the function. Common parameters include:",
    "#'   \\describe{",
    "#'     \\item{verbose}{Logical. Whether to print verbose output (default: TRUE).}",
    "#'   }",
    "#'",
    "#' @return A named list containing:",
    "#'   \\describe{",
    "#'     \\item{scRNA_data}{Modified single-cell data object with integrated screening results.}",
    "#'   }",
    "#'",
    "#' @export",
    sep = "\n",
    collapse = "\n"
  )
}

#' @keywords internal
func_body <- function(func_doc = character(), func_name = character()) {
  paste(
    sprintf(
      "%s\n%s <- function(",
      func_doc,
      func_name
    ),
    "
  matched_bulk,
  sc_data,
  phenotype,
  label_type = NULL,
  phenotype_class = c(\"binary\", \"survival\", \"continuous\"),
  ...
) {
  # Validate phenotype_class parameter
  phenotype_class <- match.arg(phenotype_class)
  
  # Extract additional arguments
  dots <- list(...)
  verbose <- dots$verbose %||% TRUE
  
  # TODO: Implement your screening logic here
  # Placeholder for modified scRNA data (replace with actual processing)
  modified_sc_data <- sc_data
  
  # Return result in expected format
  list(
    scRNA_data = modified_sc_data
  )
}",
    .sep = "\n",
    collapse = "\n"
  )
}
