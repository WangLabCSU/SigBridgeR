# ---- some scripts of functions for testing----

SaveEnv = function(
  var_list = NULL,
  path = ".",
  name = ".RData",
  preset = c("fast", "high", "balanced"),
  nthreads = min(4, parallel::detectCores())
) {
  if (length(preset) > 1 || !tolower(preset) %in% c("fast", "high", "balanced"))
    preset = "fast"
  if (is.null(var_list)) {
    all_vars = ls(envir = .GlobalEnv)
    var_list = lapply(all_vars, function(x) get(x, envir = .GlobalEnv))
    names(var_list) = all_vars
  }

  qs::qsave(
    x = var_list,
    file = file.path(path, name),
    preset = preset,
    nthreads = nthreads,
    check_hash = FALSE
  )
  cli::cli_alert_success(glue::glue(
    format(Sys.time(), "%Y/%m/%d %H:%M:%S"),
    " Environment variables saved."
  ))
  cli::cli_alert_success(glue::glue("Path: ./{name}"))
  cli::cli_alert_success(glue::glue("Working directory: {getwd()}"))
}

LoadEnv = function(path_2_file = "./.RData") {
  var_list = qs::qread(file = path_2_file)
  list2env(var_list, envir = .GlobalEnv)
  cli::cli_alert_success(sprintf(
    "%s Environment variables loaded.",
    format(Sys.time(), "%Y/%m/%d %H:%M:%S")
  ))
  cli::cli_alert_success(glue::glue("Path: {path_2_file}"))
  cli::cli_alert_success(glue::glue("Working directory: {getwd()}"))
}

LoadFunctions <- function(
  script,
  function_names,
  envir = .GlobalEnv,
  cache = TRUE
) {
  if (!exists(".script_cache", envir = .GlobalEnv)) {
    assign(".script_cache", new.env(), envir = .GlobalEnv)
  }
  cache_env <- get(".script_cache", envir = .GlobalEnv)

  script_key <- basename(script)
  if (!cache || !exists(script_key, envir = cache_env)) {
    temp_env <- new.env()
    source(script, local = temp_env)
    assign(script_key, temp_env, envir = cache_env)
  }

  loaded_funs <- lapply(
    function_names,
    function(fun_name) {
      fun_env <- get(script_key, envir = cache_env)
      if (!exists(fun_name, envir = fun_env)) {
        cli::cli_alert_danger(glue::glue(
          "Function {crayon::red(fun_name)} not found in {crayon::bold(script)}",
        ))
        stop()
      }
      get(fun_name, envir = fun_env)
    }
  )

  list2env(loaded_funs, envir = envir)
  cli::cli_alert_success(
    text = c(
      format(Sys.time(), "%Y/%m/%d %H:%M:%S"),
      glue::glue(" Function {crayon::green(function_names)} loaded"),
      glue::glue(" from {crayon::bold(script)}")
    )
  )
  return(invisible())
}

# If you see a negative optimization value, it means that is the amount of memory occupied by the function.
MemoryClean <- function(
  remove_all_objects = TRUE,
  graphics_off = FALSE,
  keep_self = TRUE,
  keep = character(),
  verbose = TRUE
) {
  mem_installed <- requireNamespace("pryr", quietly = TRUE)
  start_mem <- if (mem_installed) pryr::mem_used() else NA

  # System-level Memory Release (Linux Only)
  sys_type <- tolower(Sys.info()[["sysname"]])
  if (sys_type == "linux") {
    if (requireNamespace("Rcpp", quietly = TRUE)) {
      tryCatch(
        {
          Rcpp::sourceCpp(
            code = '
          #include <Rcpp.h>
          #include <malloc.h>
          // [[Rcpp::export]]
          void linux_malloc_trim() { malloc_trim(0); }
        '
          )
          linux_malloc_trim()
          if (verbose)
            message(
              "🔧 The Linux system-level memory release (malloc_trim) has been called."
            )
        },
        error = function(e)
          if (verbose)
            warning(
              "⚠️ System-level release failed (please install libc6-dev: sudo apt-get install libc6-dev)."
            )
      )
    } else if (verbose)
      warning(
        "⚠️ The `Rcpp` package needs to be installed to enable system-level release."
      )
  }
  # Constructing a list of preserved objects
  protected_objects <- unique(c(keep, if (keep_self) "MemoryClean"))
  # Object Cleanup
  if (remove_all_objects) {
    to_remove <- setdiff(ls(envir = .GlobalEnv), protected_objects)
    if (length(to_remove) > 0) {
      rm(list = to_remove, envir = .GlobalEnv)
      if (verbose)
        message(sprintf(
          "✅ deleted %d objects, preserved %d objects",
          length(to_remove),
          length(protected_objects)
        ))
    }
  }
  # Clean up the graphics device
  if (graphics_off) {
    grDevices::graphics.off()
    if (verbose)
      message(
        "🎨 All graphics devices have been closed."
      )
  }
  # Multiple Calls to Accelerate Release
  invisible(gc())
  invisible(gc())
  # Return Memory Report
  if (verbose && mem_installed) {
    end_mem <- pryr::mem_used()
    freed_mem <- start_mem - end_mem

    message(sprintf(
      "📊 Memory Release Report:\nBefore Release: %s MB\nAfter Release: %s MB\nReleased Amount: %s MB",
      format(start_mem / 1024^2, units = "auto"),
      format(end_mem / 1024^2, units = "auto"),
      format(freed_mem / 1024^2, units = "auto")
    ))
  }
}
