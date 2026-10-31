ValidateDirName <- function(func, ...) {
  dir_create_calls <- .vdn_find_dir_create_calls(fn_body(func))

  if (length(dir_create_calls) == 0L) {
    return(invisible())
  }

  note <- 0L
  seen <- character()

  # check the folder name of each dir.create call
  for (dir_call in dir_create_calls) {
    path_arg <- .vdn_path_arg(dir_call)
    folder_name <- if (is.null(path_arg)) {
      NULL
    } else {
      .vdn_folder_name(path_arg, func)
    }

    # for complex expressions, record a warning and continue
    if (is.null(folder_name)) {
      cli::cli_warn("Failed to get folder name for dir.create call")
      next
    }

    # * the same folder is only reported once, e.g. when both branches of an
    # * `if` expression create it
    if (folder_name %in% seen) {
      next
    }
    seen <- c(seen, folder_name)

    # "_res"
    if (.vdn_has_res_suffix(folder_name)) {
      validate_success(
        "Folder name '{folder_name}' ends with '_res'",
        .envir = current_env()
      )

      next
    }

    # suggested folder name
    suggested_folder_name <- paste0(folder_name, "_res")

    validate_note("Folder name does not end with '_res'")
    validate_explain(
      c(
        "{cli::symbol$bullet} Current folder name: '{folder_name}'",
        "{cli::symbol$bullet} Recommended: '{suggested_folder_name}'",
        "Result ('_res') folders stay apart from source code"
      ),
      .envir = current_env()
    )
    note <- note + 1L
  }

  if (note > 0L) {
    return(list(note = note))
  }

  invisible()
}

# Recursively collects every `dir.create()` call in `expr`.
# `is_call()` also accepts formulas, which are rejected by `call_name()` and
# `call_args()`, hence the `safe_*()` helpers.
.vdn_find_dir_create_calls <- function(expr) {
  calls <- list()

  if (is_call(expr)) {
    if (identical(safe_call_name(expr), "dir.create")) {
      calls <- list(expr)
    }
    # recursively traverse sub-expressions
    for (arg in safe_call_args(expr)) {
      calls <- c(calls, .vdn_find_dir_create_calls(arg))
    }
  } else if (is_pairlist(expr)) {
    for (elt in expr) {
      calls <- c(calls, .vdn_find_dir_create_calls(elt))
    }
  }

  calls
}

# Returns the expression matched to the `path` argument of a `dir.create()`
# call. `path` is the first formal of `dir.create()`, so the first untagged
# argument is matched to it positionally.
.vdn_path_arg <- function(dir_call) {
  args <- safe_call_args(dir_call)

  if (length(args) == 0L) {
    return(NULL)
  }

  arg_names <- names(args)

  if (!is.null(arg_names) && "path" %in% arg_names) {
    return(args[["path"]])
  }

  untagged <- if (is.null(arg_names)) {
    seq_along(args)
  } else {
    which(arg_names == "")
  }

  if (length(untagged) == 0L) {
    return(NULL)
  }

  args[[untagged[1L]]]
}

# Resolves the folder name created by a `dir.create()` call.
# @param path_arg The expression matched to the `path` argument.
# @param func The function being validated, used to evaluate variables.
# @return A length-1 character, or `NULL` when the name cannot be determined.
.vdn_folder_name <- function(path_arg, func) {
  folder_name <- if (is_syntactic_literal(path_arg)) {
    path_arg
  } else if (is_symbol(path_arg)) {
    # the folder name may be stored in a variable of the function
    tryCatch(
      get_var_value(var_name = as.character(path_arg), func = func),
      error = function(e) NULL
    )
  } else {
    return(NULL)
  }

  if (
    !is.atomic(folder_name) ||
      length(folder_name) != 1L ||
      is.na(folder_name)
  ) {
    return(NULL)
  }

  # `dir.create("foo_res/")` still names a folder suffixed with `_res`
  gsub("[/\\\\]+$", "", as.character(folder_name))
}

# Does the folder name end with the `_res` suffix?
.vdn_has_res_suffix <- function(folder_name) {
  grepl("_res$", folder_name)
}
