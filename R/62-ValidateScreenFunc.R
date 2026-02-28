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

  if (rlang::is_installed("tictoc")) {
    tictoc::tic()
  }

  cli::cli_text(sprintf(
    "Start at %s",
    format(Sys.time(), "%Y/%m/%d %H:%M:%S")
  ))
  cli::cli_text("\n")

  inputs_check <- ValidateArgsInputs(func = func, ...)
  cli::cli_text("\n")

  verbose_check <- ValidateArgsVerbose(func = func, ...)
  cli::cli_text("\n")

  syntax_check <- ValidateArgsSyntax(func = func, ...)
  cli::cli_text("\n")

  return_check <- ValidateReturn(func = func, ...)
  cli::cli_text("\n")

  dir_name_check <- ValidateDirName(func = func, ...)
  cli::cli_text("\n")

  if (rlang::is_installed("tictoc")) {
    end_time <- tictoc::toc(quiet = TRUE)
    cli::cli_text(tictoc::toc.outmsg(end_time$tic, end_time$toc, "Duration"))
    cli::cli_text("\n")
  }

  ValidateBadge(inputs_check, verbose_check, syntax_check, return_check)
  cli::cli_text("\n")

  invisible(TRUE)
}


#' @keywords internal
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
  has_dots <- "..." %in% func_args

  # Check required parameters
  missing_args <- setdiff(required_args, func_args)
  if (length(missing_args) == 0) {
    validate_success("All input arguments explicitly specified")
    return(list(error = 0))
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
    theme = list(div = list("margin-left" = 2))
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

  list(error = 1)
}

#' @keywords internal
ValidateArgsVerbose <- function(func, ...) {
  func_formals <- formals(func)
  func_args <- names(func_formals)

  verbose_used <- any(grepl("verbose", as.character(rlang::fn_body(func))))
  has_dots <- "..." %in% func_args
  has_verbose <- "verbose" %in% func_args || has_dots

  if (!has_verbose) {
    validate_note("Verbose control not supported")
    validate_explain(
      "Consider adding `verbose` control to ease error tracing"
    )
    return(list(note = 1))
  } else if (!verbose_used) {
    validate_warn("Verbose control supported but not used")
    if (has_dots) {
      validate_explain(
        "{cli::symbol$bullet} When using `...` to accept extra arguments, also support `verbose`"
      )
    }
    return(list(warn = 1))
  }

  validate_success("Verbose control supported")

  invisible()
}

#' @keywords internal
ValidateArgsSyntax <- function(func, ...) {
  if (rlang::is_installed(c("tidycheckUsage", "knitr"))) {
    report <- tidycheckUsage::tidycheckUsage(func)
    if (!is.null(report)) {
      validate_error("Syntax error in function")
      report$file <- NULL
      report$path <- NULL
      report$fun <- NULL
      print(knitr::kable(
        report,
        format = "pipe",
        align = "c",
        row.names = FALSE
      ))
      return(list(error = 1))
    }
    validate_success("Syntax check passed")
    return(invisible())
  }

  if (rlang::check_installed("codetools")) {
    code_report <- utils::capture.output(
      codetools::checkUsage(
        func,
        name = "Syntax",
        report = validate_explain
      ),
      type = "message"
    )
    if (!is.null(code_report)) {
      cli::cli_par()
      validate_error("Syntax error in function")
      for (msg in code_report) {
        validate_explain(msg)
      }
      cli::cli_end()

      error <- 1
    } else {
      validate_success("Syntax check passed")
      error <- 0
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
    if (length(suspicious_globals) != 0) {
      validate_warn("Suspicious global variables")
      validate_explain("Undefined variables detected in function:")

      cli::cli_div(
        theme = list(div = list("margin-left" = 4))
      )
      cli::cli_text("{.val {suspicious_globals}}")
      cli::cli_end()
      warn <- 1
    } else {
      validate_success("No suspicious global variables")
      warn <- 0
    }

    return(list(error = error, warn = warn))
  }

  cli::cli_warn(
    "Syntax check not supported due to missing dependencies:\\
    \n  {.pkg tidycheckUsage & knitr} or {.pkg codetools}"
  )
}

#' @keywords internal
ValidateReturn <- function(func, ...) {
  func_body <- rlang::fn_body(func)

  find_last_return <- function(expr) {
    if (!is.call(expr)) {
      expr
    }

    if (as.character(expr[[1]]) == "return") {
      return(expr[[2]])
    }

    if (as.character(expr[[1]]) == "{") {
      last_expr <- expr[[length(expr)]]
      return(find_last_return(last_expr))
    }

    if (as.character(expr[[1]]) == "if") {
      if (length(expr) >= 3) {
        return(find_last_return(expr[[3]]))
      }
    }

    expr
  }

  last_return <- find_last_return(func_body)

  is_list_creation <- FALSE
  has_scRNA_data <- FALSE

  if (is.call(last_return)) {
    func_name <- as.character(last_return[[1]])

    if (func_name == "list") {
      is_list_creation <- TRUE

      arg_names <- names(last_return)[-1]

      if (is.null(arg_names) || all(arg_names == "")) {
        for (i in 2:length(last_return)) {
          arg <- last_return[[i]]
          if (is.call(arg) && as.character(arg[[1]]) == "=") {
            arg_name <- as.character(arg[[2]])
            if (arg_name == "scRNA_data") {
              has_scRNA_data <- TRUE
              break
            }
          }
        }
      } else {
        has_scRNA_data <- "scRNA_data" %in% arg_names
      }
    }
  }

  if (!is_list_creation) {
    validate_error("Return value is not a list")
    validate_explain_speaker(
      message = c(
        "recommended to be the first element of the return value",
        "Should be of class {.cls Seurat}"
      ),
      message_name = "scRNA_data",
      type = "error"
    )
    return(list(error = 1))
  } else if (!has_scRNA_data) {
    validate_error("Return value does not have `scRNA_data` slot")
    validate_explain_speaker(
      message = "Should be of class {.cls Seurat}",
      message_name = "scRNA_data",
      type = "error"
    )
    return(list(error = 1))
  } else {
    validate_success("Return value is a list with `scRNA_data` slot")
  }

  invisible()
}

#' @keywords internal
ValidateDirName <- function(func, ...) {
  func_body <- rlang::fn_body(func)

  # 递归查找所有 dir.create 调用
  find_dir_create_calls <- function(expr) {
    calls <- list()
    if (rlang::is_call(expr)) {
      if (rlang::call_name(expr) == "dir.create") {
        calls <- list(expr)
      }
      # 递归遍历子表达式
      for (arg in rlang::call_args(expr)) {
        calls <- c(calls, find_dir_create_calls(arg))
      }
    } else if (rlang::is_pairlist(expr)) {
      for (elt in expr) {
        calls <- c(calls, find_dir_create_calls(elt))
      }
    }
    calls
  }

  dir_create_calls <- find_dir_create_calls(func_body)

  if (length(dir_create_calls) == 0) {
    return(invisible())
  }

  note <- 0L

  # 检查每个 dir.create 调用的文件夹名
  for (call in dir_create_calls) {
    folder_arg <- rlang::call_args(call)[[1]]

    # 尝试获取文件夹名
    if (rlang::is_syntactic_literal(folder_arg)) {
      folder_name <- as.character(folder_arg)
    } else if (rlang::is_symbol(folder_arg)) {
      # 如果是符号，尝试获取变量名
      folder_name <- as.character(folder_arg)
    } else {
      # 对于复杂表达式，记录警告但继续
      validate_note(
        "Cannot statically analyze folder name from complex expression"
      )
      next
    }

    # "_res"
    if (!grepl("_res", folder_name)) {
      validate_note("Folder name does not end with '_res'")
      validate_explain(
        c(
          glue::glue(
            "{cli::symbol$bullet} Current folder name: '{folder_name}'"
          ),
          glue::glue(
            "{cli::symbol$bullet} Recommended: Use '_res' suffix (e.g., '{folder_name}_res') to clearly distinguish the result (intermediate) folder from the source code"
          )
        )
      )
      note <- note + 1L
    } else {
      validate_success(
        "Folder name '{folder_name}' ends with '_res'"
      )
    }
  }

  if (note > 0) {
    return(list(note = note))
  }
  invisible()
}


#' @keywords internal
ValidateBadge <- function(...) {
  merge_check_res <- function(...) {
    lists <- rlang::list2(...)

    all_names <- c("error", "warn", "note")
    lists <- lapply(lists, function(x) {
      missing <- setdiff(all_names, names(x))
      if (length(missing)) {
        x[missing] <- 0
      }
      x[all_names]
    })

    purrr::reduce(lists, ~ purrr::map2(.x, .y, `+`))
  }
  check_res <- merge_check_res(...)

  symbol <- function(value = integer(), name = character()) {
    t <- cli::symbol$tick
    x <- cli::symbol$cross

    if (value == 0) {
      cli::cli_fmt(cli::cli_text(cli::col_green("{value} {name} {t}")))
    } else if (name == "note") {
      if (value == 1) {
        cli::cli_fmt(cli::cli_text(cli::col_blue("{value} {name} {x}")))
      } else {
        cli::cli_fmt(cli::cli_text(cli::col_blue("{value} {name}s {x}")))
      }
    } else if (value == 1) {
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
