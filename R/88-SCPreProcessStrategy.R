#' @title Preprocessing Strategy Registry for Single-Cell Workflows
#'
#' @description
#' An environment that maps single-letter codes to standardized Seurat preprocessing
#' operations. Each entry is a function of the form \code{function(object, params)}
#' that wraps a core Seurat step (e.g., normalization, PCA, clustering) using
#' \code{rlang::exec()} for safe parameter injection.
#'
#' This registry enables compact, composable pipeline definitions (e.g., via character
#' strings like \code{"onvps"} for "Create → Normalize → VariableFeatures → Scale → PCA")
#' and supports both standard (\code{NormalizeData}) and SCTransform-based workflows.
#'
#' @section Available Operations:
#' \describe{
#'   \item{\code{o}}{Create Seurat object from count matrix (\code{SeuratObject::CreateSeuratObject})}
#'   \item{\code{n}}{Normalize data (\code{Seurat::NormalizeData})}
#'   \item{\code{v}}{Find variable features (\code{Seurat::FindVariableFeatures})}
#'   \item{\code{s}}{Scale data (\code{Seurat::ScaleData})}
#'   \item{\code{p}}{Run PCA (\code{Seurat::RunPCA})}
#'   \item{\code{c}}{Find clusters (\code{Seurat::FindClusters})}
#'   \item{\code{e}}{Find neighbors (\code{Seurat::FindNeighbors})}
#'   \item{\code{t}}{Run t-SNE (\code{Seurat::RunTSNE})}
#'   \item{\code{u}}{Run UMAP (\code{Seurat::RunUMAP})}
#'   \item{\code{r}}{SCTransform (\code{Seurat::SCTransform})}
#'   \item{\code{i}}{Integrate layers (\code{Seurat::IntegrateLayers})}
#' }
#'
#' @section Usage Notes:
#' \itemize{
#'   \item These functions are **not intended for direct interactive use**. They are
#'         internal building blocks for workflow engines (i.e., \code{SCPreProcess}).
#'   \item Each function expects \code{object} (a matrix or Seurat object) and
#'         \code{params} (a named list of arguments passed via \code{!!!} splicing).
#'   \item You can access any operation via \code{SCPreProcessStrategy$letter},
#'         but doing so bypasses pipeline validation and error handling.
#'   \item To add more operations, use \code{RegisterSeuratMethod()}
#' }
#'
#' @family single_cell_preprocess
#' @export
SCPreProcessStrategy <- rlang::new_environment(
  list(
    o = function(object, params) {
      rlang::exec(
        SeuratObject::CreateSeuratObject,
        counts = object,
        !!!params
      )
    },
    n = function(object, params) {
      rlang::exec(
        Seurat::NormalizeData,
        object = object,
        !!!params
      )
    },
    v = function(object, params) {
      rlang::exec(
        Seurat::FindVariableFeatures,
        object = object,
        !!!params
      )
    },
    s = function(object, params) {
      rlang::exec(
        Seurat::ScaleData,
        object = object,
        !!!params
      )
    },
    p = function(object, params) {
      rlang::exec(
        Seurat::RunPCA,
        object = object,
        !!!params
      )
    },
    c = function(object, params) {
      rlang::exec(
        Seurat::FindClusters,
        object = object,
        !!!params
      )
    },
    e = function(object, params) {
      rlang::exec(
        Seurat::FindNeighbors,
        object = object,
        !!!params
      )
    },
    t = function(object, params) {
      rlang::exec(
        Seurat::RunTSNE,
        object = object,
        !!!params
      )
    },
    u = function(object, params) {
      rlang::exec(
        Seurat::RunUMAP,
        object = object,
        !!!params
      )
    },
    r = function(object, params) {
      rlang::exec(
        Seurat::SCTransform,
        object = object,
        !!!params
      )
    },
    i = function(object, params) {
      rlang::exec(
        Seurat::IntegrateLayers,
        object = object,
        !!!params
      )
    }
  )
)
