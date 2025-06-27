# --------------------- Log Output Function ---------------------

#' @description This function is used for a quick check of the execution status. If more detailed information is needed, saving the R script console information would be a better choice.
#' @param ... The content of the log message, here auto-using `glue()` to generate message.
#' @param log_tpye The type of log message. It can be one of "success", "danger", "info", "warn".
#' @param show_log Whether to show log message in console.
#' @param log_operation The operation to be performed on log content. It can be one of "reset", "to_file".
#' @examples
#' # initialize log
#' LogOutput("This is a log message")
#' # reset log
#' LogOutput(log_operation = "reset")
#' # add complex log with a log type
#' LogOutput("Processed {n} items", n = 100, log_type = "success")
#' # save log to file
#' LogOutput(log_operation = ToFile(file="./mylog.log"))
LogOutput <- local({
  # Private environment to store log content and configurations
  .log_env <- new.env(parent = emptyenv())
  .log_env$content <- character(0)
  .log_env$cli_available <- requireNamespace("cli", quietly = TRUE)

  TimeStamp <- function() format(Sys.time(), "%Y/%m/%d %H:%M:%S")

  GenerateMessage <- function(..., env = parent.frame()) {
    if (length(list(...)) > 0) {
      msg <- glue::glue(..., .envir = env)
      glue::glue("[{TimeStamp()}] {trimws(msg)}")
    } else {
      glue::glue("[{TimeStamp()}] Log initialized")
    }
  }

  HandleOperation <- function(operation) {
    if (inherits(operation, "to_file")) {
      tryCatch(
        {
          writeLines(.log_env$content, operation$filename)
          success_msg <- glue::glue("Log saved to {operation$filename}")
          .log_env$content <<- c(.log_env$content, success_msg)
        },
        error = function(e) warning("File write failed: ", e$message)
      )
    } else if (identical(tolower(operation), "reset")) {
      .log_env$content <<- character(0)
    } else {
      warning("Unsupported operation type: ", toString(class(operation)))
    }
  }

  function(
    ...,
    log_type = NULL,
    show_log = TRUE,
    log_operation = NULL
  ) {
    if (!is.null(log_operation)) {
      HandleOperation(log_operation)
      return(invisible(.log_env$content))
    }

    new_msg <- GenerateMessage(..., env = parent.frame())

    if (!is.null(log_type) && .log_env$cli_available) {
      new_msg <- switch(
        EXPR = log_type,
        success = glue::glue("✅", new_msg),
        danger = glue::glue("❌", new_msg),
        info = glue::glue("ℹ️", new_msg),
        warn = glue::glue("⚠️", new_msg),
        inspire = glue::glue("✨", new_msg),
        setting = glue::glue("🔧", new_msg),
        start = glue::glue("▶️", new_msg),
        end = glue::glue("⏹️", new_msg),
        new_msg
      )
    }

    .log_env$content <<- c(.log_env$content, new_msg)
    # print log
    if (show_log) {
      if (.log_env$cli_available && !is.null(log_type)) {
        cli::cli_alert(new_msg)
      } else {
        message(new_msg)
      }
    }

    return(invisible(.log_env$content))
  }
})

ToFile <- function(file = "mylog.log") {
  # .txt also available
  structure(list(filename = file), class = "to_file")
}
