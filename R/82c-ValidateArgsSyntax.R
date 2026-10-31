ValidateArgsSyntax <- function(func, ...) {
  if (is_installed("codetools")) {
    code_report <- utils::capture.output(
      codetools::checkUsage(
        func,
        name = "Syntax",
        report = validate_explain,
        all = TRUE
      ),
      type = "message"
    )
    if (length(code_report) > 0L) {
      cli::cli_par()
      validate_error("Syntax error in function")
      for (msg in code_report) {
        validate_explain(msg)
      }
      cli::cli_end()

      error <- 1L
    } else {
      validate_success("Syntax check passed")
      error <- 0L
    }

    globals <- codetools::findGlobals(func)

    safe_symbols <- unique(c(
      ls(envir = baseenv()),
      ls(envir = as.environment("package:stats")),
      ls(envir = as.environment("package:utils")),
      ls(envir = as.environment("package:methods")),
      ".Primitive",
      ".Call",
      ".C",
      ".Fortran",
      ".External"
    ))
    suspicious_globals <- setdiff(globals, safe_symbols)
    if (length(suspicious_globals) != 0L) {
      validate_warn("Suspicious global variables")
      validate_explain("Undefined variables detected in function:")

      cli::cli_div(
        theme = list(div = list("margin-left" = 4L))
      )
      cli::cli_text("{.val {suspicious_globals}}")
      cli::cli_end()
      warn <- 1L
    } else {
      validate_success("No suspicious global variables")
      warn <- 0L
    }

    return(list(error = error, warn = warn))
  }

  cli::cli_warn(
    "Syntax check not supported due to missing dependencies: \
      {.pkg codetools}"
  )

  # * keep the result shape stable for `ValidateBadge()`
  list(error = 0L, warn = 0L)
}
