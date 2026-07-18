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

  invisible()
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# suppress R CMD CHECK NOTE about global variables when using tidyverse or data.table
utils::globalVariables(c(
  "PC",
  "PC1",
  "PC2",
  "sd",
  "base_key",
  ".",
  "suffix",
  "max_suffix",
  "row_id",
  "value",
  "col_name",
  "condition",
  "batch",
  "composite_score",
  "variance_stability",
  "marker_signal",
  "dropout_robustness",
  "method",
  "DEGAS.model_type",
  "DEGAS.architecture",
  "DEGAS.ff_depth",
  "DEGAS.bag_depth",
  "DEGAS.seed",
  "Feature",
  "..duplicate_cols",
  "..vote_cols",
  "n",
  "Total",
  "Fraction",
  "sets",
  "count",
  "Variance",
  "Rank_Label",
  "Rank_Score",
  "Reject",
  "TiRank"
))
