# nocov start

# ? Package startup messages
.onAttach <- function(libname, pkgname) {
  if (startup_message_allowed()) {
    pkg_version <- utils::packageVersion(pkgname)
    end <- Sys.time()
    elapsed <- cli_timestamp_formatter(as.numeric(end - .time_record$start))
    cli::cli_alert_success(sprintf(
      "{.pkg {pkgname}} v{pkg_version} loaded {cli::col_grey('[%s]')}",
      elapsed
    ))
  }
}

.onLoad <- function(libname, pkgname) {
  if (startup_message_allowed()) {
    .time_record$start <- Sys.time()
  }

  # * default options
  op <- options()
  op_pkg <- list(
    SigBridgeR.verbose = TRUE,
    SigBridgeR.seed = 123L,
    SigBridgeR.timeout = 180L,
    SigBridgeR.nthreads = 2L,
    IDConverter.datapath = system.file("extdata", package = "IDConverter")
  )

  toset <- !(names(op_pkg) %chin% names(op))
  if (any(toset)) {
    options(op_pkg[toset])
  }

  # * Python module
  reticulate::py_require(
    packages = c(
      "celltypist",
      "anndata",
      "pandas"
    ),
    python_version = ">=3.6"
  )
  python_path <- system.file(
    "python",
    package = pkgname
  )
  PyModule$celltypist <- reticulate::import_from_path(
    module = "celltypist_annotate",
    path = python_path,
    convert = FALSE,
    delay_load = TRUE
  )

  # * S7 OOP system
  S7::methods_register()
  invisible()
}

.time_record <- new.env()

cli_timestamp_formatter <- function(x) {
  stopifnot(
    length(x) == 1L,
    is.numeric(x),
    !is.na(x),
    x >= 0L
  )

  if (x < 1L) {
    sprintf("%.2f ms", x * 1000L)
  } else if (x < 60L) {
    sprintf("%.2f s", x)
  } else if (x < 3600L) {
    sprintf("%d min %.2f s", floor(x / 60L), x %% 60L)
  } else {
    sprintf(
      "%d h %02d min %.2f s",
      floor(x / 3600L),
      floor((x %% 3600L) / 60L),
      x %% 60L
    )
  }
}

`%||%` <- function(x, y) if (is.null(x)) y else x

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
