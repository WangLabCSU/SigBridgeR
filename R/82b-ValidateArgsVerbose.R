ValidateArgsVerbose <- function(func, ...) {
  func_formals <- formals(func)
  func_args <- names(func_formals)

  verbose_used <- any(grepl("verbose", as.character(fn_body(func))))
  has_dots <- "..." %chin% func_args
  has_verbose <- "verbose" %chin% func_args || has_dots

  if (!has_verbose) {
    validate_note("Verbose control not supported")
    validate_explain(
      "Consider adding `verbose` control to ease error tracing"
    )
    return(list(note = 1L))
  } else if (!verbose_used) {
    validate_warn("Verbose control supported but not used")
    if (has_dots) {
      validate_explain(
        "{cli::symbol$bullet} When using `...` to accept extra arguments, also support `verbose`"
      )
    }
    return(list(warn = 1L))
  }

  validate_success("Verbose control supported")

  invisible()
}
