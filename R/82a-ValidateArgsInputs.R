ValidateArgsInputs <- function(func, ...) {
  required_args <- c(
    "sc_data",
    "matched_bulk",
    "phenotype",
    "label_type",
    "phenotype_class"
  )
  func_formals <- formals(func)
  func_args <- names(func_formals)
  has_dots <- "..." %chin% func_args

  # Check required parameters
  missing_args <- setdiff(required_args, func_args)
  if (length(missing_args) == 0L) {
    validate_success("All input arguments explicitly specified")
    return(list(error = 0L))
  }

  validate_error("Missing required arguments")
  validate_explain("More arguments should be added:")
  reason <- list(
    sc_data = "A fully preprocessed Seurat object",
    matched_bulk = c(
      "Bulk RNA-seq matrix (gene * samples):",
      "Samples match `phenotype` in number/order;",
      "Senes overlap with `sc_data`"
    ),
    phenotype = c(
      "Phenotype: a named vector or data.frame, names/rownames match `matched_bulk`",
      "For binary/continuous: named vector recommended",
      "For survival: data.frame with `time` (1st col) and `status` (2nd col) recommended."
    ),
    label_type = "Labeling phenotype-associated cell with real study identifiers",
    phenotype_class = c(
      "Phenotype types:",
      "Must be one or more of `binary`, `continuous` and `survival`"
    )
  )
  cli::cli_div(
    id = "missing_args",
    theme = list(div = list("margin-left" = 2L))
  )
  purrr::walk(
    missing_args,
    ~ validate_explain_speaker(
      reason[[.x]],
      .x,
      "error"
    )
  )
  cli::cli_end(id = "missing_args")

  list(error = 1L)
}
