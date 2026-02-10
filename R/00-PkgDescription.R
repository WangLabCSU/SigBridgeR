# ? Package Description

#' @title SigBridgeR: Integrative Framework and Toolkit for Single-Cell Screening of Phenotype-Associated Cells
#'
#' @description
#' SigBridgeR is an integrative toolkit designed to identify phenotype-associated cell subpopulations by combining phenotype(e.g. survival, drug sensitivity), bulk expression and single-cell RNA-seq data. It leverages multiple algorithms to robustly link cell features with clinical or functional phenotypes. The package provides a unified pipeline for cross-modal data analysis, enabling the discovery of biologically and clinically relevant cell states in heterogeneous samples.
#'
#' @section Main functions:
#' The package includes multiple algorithms for integrative analysis of single-cell and bulk data to identify phenotype-associated cell populations. (see function \code{\link{Screen}})
#'
#' @section Data requirements (their pre-processing is also provided in the package):
#' - Single-cell RNA-seq data (Seurat object format)
#' - Bulk expression data
#' - Phenotype data (survival, drug sensitivity, etc.)
#'
#' @section Tutorial:
#' For a detailed tutorial, please visit: \url{https://wanglabcsu.github.io/SigBridgeR/}
#'
#' @section License:
#' SigBridgeR is licensed under the GPL version 3.0. Please see the LICENSE file for details.
#'
#'
#' @author Yuxi Yang \email{15364051195@163.com} ORCID: 0009-0006-1329-1224 (creator, author)
#'
#'
#' @docType package
#' @name SigBridgeR-package
#' @aliases SigBridgeR
#' @keywords internal
#'
"_PACKAGE"
