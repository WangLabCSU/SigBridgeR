# ? Package startup messages
.onAttach <- function(libname, pkgname) {
  pkg_version <- utils::packageVersion(pkgname)

  msg <- cli::cli_fmt(cli::cli_alert_success(
    "{.pkg {pkgname}} v{pkg_version} loaded"
  ))
  packageStartupMessage(msg)
  invisible()
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
