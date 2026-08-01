# nocov start

# ? Package startup messages
.onAttach <- function(libname, pkgname) {
  pkg_version <- utils::packageVersion(pkgname)

  startup_spinner(
    expr = invisible(),
    pkgname = pkgname,
    pkg_version = pkg_version
  )
}

.onLoad <- function(libname, pkgname) {
  # default options
  op <- options()
  op_pkg <- list(
    SigBridgeR.verbose = TRUE,
    SigBridgeR.seed = 123L,
    SigBridgeR.timeout = 180L,
    IDConverter.datapath = system.file("extdata", package = "IDConverter")
  )

  toset <- !(names(op_pkg) %chin% names(op))
  if (any(toset)) {
    options(op_pkg[toset])
  }
  S7::methods_register()
  invisible()
}

`%||%` <- function(x, y) if (is.null(x)) y else x

startup_spinner <- function(expr, pkgname, pkg_version) {
  if (!is_interactive() || !startup_message_allowed()) {
    return(force(expr))
  }

  id <- cli::cli_progress_step(
    msg = "{.pkg {pkgname}} v{pkg_version} loading",
    msg_done = "{.pkg {pkgname}} v{pkg_version} loaded",
    msg_failed = "{.pkg {pkgname}} v{pkg_version} fail to load",
    spinner = TRUE
  )

  on.exit(
    cli::cli_progress_done(id = id),
    add = TRUE
  )

  force(expr)
}

startup_message_allowed <- function() {
  allowed <- FALSE

  withRestarts(
    {
      signalCondition(structure(
        list(message = ".__startup_probe__."),
        class = c(
          "packageStartupProbe",
          "packageStartupMessage",
          "condition"
        )
      ))
      allowed <- TRUE
    },
    muffleMessage = function() NULL
  )

  allowed
}
