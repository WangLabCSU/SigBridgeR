#' @title Preprocessing Strategy Registry for Single-Cell Workflows
#'
#' @description
#' An environment that maps single-letter codes to standardized Seurat preprocessing
#' operations. Each entry is a function of the form \code{function(...)}
#' that wraps a core Seurat step (e.g., normalization, PCA, clustering)
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
#'   \item You can access any operation via \code{SCPreProcessStrategy$letter},
#'         but doing so bypasses pipeline validation and error handling.
#'   \item To add more operations, use \code{RegisterSeuratMethod()}
#' }
#'
#' @family single_cell_preprocess
#' @export
SCPreProcessStrategy <- rlang::new_environment(
  list(
    o = function(...) {
      SeuratObject::CreateSeuratObject(...)
    },
    n = function(...) {
      Seurat::NormalizeData(...)
    },
    v = function(...) {
      Seurat::FindVariableFeatures(...)
    },
    s = function(...) {
      Seurat::ScaleData(...)
    },
    p = function(...) {
      dots <- rlang::list2(...)
      if (dots$features == "all") {
        dots$features <- rownames(dots[[1]])
      }
      rlang::exec(Seurat::RunPCA, !!!dots)
    },
    c = function(...) {
      Seurat::FindClusters(...)
    },
    e = function(...) {
      Seurat::FindNeighbors(...)
    },
    t = function(...) {
      Seurat::RunTSNE(...)
    },
    u = function(...) {
      Seurat::RunUMAP(...)
    },
    r = function(...) {
      Seurat::SCTransform(...)
    },
    i = function(...) {
      Seurat::IntegrateLayers(...)
    }
  )
)
