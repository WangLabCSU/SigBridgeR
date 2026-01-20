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
    SigBridgeR.timeout = 180L
  )

  toset <- !(names(op_pkg) %chin% names(op))
  if (any(toset)) {
    options(op_pkg[toset])
  }

  # register screening method
  RegisterScreenMethod(
    "Scissor" = DoScissor,
    "scPAS" = DoscPAS,
    supported_phenotypes = c("binary", "survival", "continuous"),
    parameter_mapper = function(params) {
      params$family <- switch(
        params$phenotype_class,
        "binary" = "binomial",
        "survival" = "cox",
        "continuous" = "gaussian"
      )
      params
    },
    registry = ScreenStrategy,
    verbose = FALSE
  )

  RegisterScreenMethod(
    "scPP" = DoscPP,
    "DEGAS" = DoDEGAS,
    supported_phenotypes = c("binary", "survival", "continuous"),
    parameter_mapper = NULL,
    registry = ScreenStrategy,
    verbose = FALSE
  )

  RegisterScreenMethod(
    "scAB" = DoscAB,
    supported_phenotypes = c("binary", "survival"),
    parameter_mapper = NULL,
    registry = ScreenStrategy,
    verbose = FALSE
  )

  RegisterScreenMethod(
    "LP_SGL" = DoLP_SGL,
    supported_phenotypes = c("binary", "survival", "continuous"),
    parameter_mapper = function(params) {
      params$family <- switch(
        params$phenotype_class,
        "binary" = 'logit',
        "survival" = 'cox',
        "continuous" = 'linear'
      )
      params
    },
    registry = ScreenStrategy,
    verbose = FALSE
  )

  RegisterScreenMethod(
    "PIPET" = DoPIPET,
    supported_phenotypes = c("binary", "continuous", "survival"),
    parameter_mapper = NULL,
    registry = ScreenStrategy,
    verbose = FALSE
  )

  invisible()
}
