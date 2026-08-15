# ? Package Description
# nocov start

#' @title SigBridgeR: Integrative Framework and Toolkit for Single-Cell Screening of Phenotype-Associated Cells
#'
#' @description
#' SigBridgeR is an integrative toolkit designed to identify phenotype-associated cell subpopulations by combining phenotype(e.g. survival, drug sensitivity), bulk expression and single-cell RNA-seq data. It leverages multiple algorithms to robustly link cell features with clinical or functional phenotypes. The package provides a unified pipeline for cross-modal data analysis, enabling the discovery of biologically and clinically relevant cell states in heterogeneous samples.
#'
#' @section Main functions:
#' The package includes multiple algorithms for integrative analysis of single-cell and bulk data to identify phenotype-associated cell populations. (see function \code{\link{Screen}})
#'
#' @section Data requirements:
#' - Single-cell RNA-seq data (Seurat object format)
#' - Bulk expression data
#' - Phenotype data (survival, drug sensitivity, etc.)
#'
#' @section Tutorial:
#' For a detailed tutorial, please visit: \url{https://wanglabcsu.github.io/SigBridgeR/}
#'
#' @section Citation:
#' Use `citation("SigBridgeR")` to cite SigBridgeR in publications.
#'
#' @section License:
#' SigBridgeR is licensed under the GPL version 3.0. Please see the LICENSE file for details.
#'
#'
#' @author `r rd_authors()`
#'
#' @docType package
#' @name SigBridgeR-package
#' @aliases SigBridgeR
#' @keywords internal
#'
"_PACKAGE"

rd_authors <- function() {
  # Read Authors@R from DESCRIPTION dynamically using base R
  # roxygen2 runs in the package root, so DESCRIPTION is in the working dir
  dcf <- read.dcf("DESCRIPTION")
  authors <- eval(parse(text = dcf[1, "Authors@R"]))

  # Format each author as a roxygen @author entry
  lines <- vapply(
    authors,
    function(p) {
      name <- paste(p$given, p$family)
      email <- p$email
      orcid <- p$comment["ORCID"]

      if (length(email) && nzchar(email)) {
        name <- sprintf("%s \\email{%s}", name, email)
      }
      if (length(orcid) && nzchar(orcid)) {
        name <- sprintf("%s (ORCID: %s)", name, orcid)
      }
      name
    },
    character(1)
  )

  paste(lines, collapse = "\n")
}
