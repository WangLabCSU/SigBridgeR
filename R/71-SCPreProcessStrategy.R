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
#'
#' * **`o`**: Create Seurat object from count matrix (`SeuratObject::CreateSeuratObject`)
#' * **`n`**: Normalize data (`Seurat::NormalizeData`)
#' * **`v`**: Find variable features (`Seurat::FindVariableFeatures`)
#' * **`s`**: Scale data (`Seurat::ScaleData`)
#' * **`p`**: Run PCA (`Seurat::RunPCA`)
#' * **`c`**: Find clusters (`Seurat::FindClusters`)
#' * **`e`**: Find neighbors (`Seurat::FindNeighbors`)
#' * **`t`**: Run t-SNE (`Seurat::RunTSNE`)
#' * **`u`**: Run UMAP (`Seurat::RunUMAP`)
#' * **`r`**: SCTransform (`Seurat::SCTransform`)
#' * **`i`**: Integrate layers (`Seurat::IntegrateLayers`)
#'
#' @section Usage Notes:
#'
#' * These functions are **not intended for direct interactive use**. They are
#'   internal building blocks for workflow engines (i.e., `SCPreProcess`).
#' * You can access any operation via `SCPreProcessStrategy$letter`,
#'   but doing so bypasses pipeline validation and error handling.
#' * To add more operations, use `RegisterSeuratMethod()`
#'
#' @family single_cell_preprocess
#' @export
SCPreProcessStrategy <- new_environment(
  list(
    o = SeuratObject::CreateSeuratObject,
    n = Seurat::NormalizeData,
    v = Seurat::FindVariableFeatures,
    s = Seurat::ScaleData,
    p = function(...) {
      dots <- rlang::list2(...)
      if (is.null(dots$features)) {
        return(Seurat::RunPCA(...))
      } else if (dots$features == "all") {
        dots$features <- rownames(dots[[1]])
      }
      do.call(Seurat::RunPCA, dots)
    },
    c = Seurat::FindClusters,
    e = Seurat::FindNeighbors,
    t = Seurat::RunTSNE,
    u = Seurat::RunUMAP,
    r = Seurat::SCTransform,
    i = Seurat::IntegrateLayers
  )
)
