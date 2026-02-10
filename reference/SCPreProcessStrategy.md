# Preprocessing Strategy Registry for Single-Cell Workflows

An environment that maps single-letter codes to standardized Seurat
preprocessing operations. Each entry is a function of the form
`function(...)` that wraps a core Seurat step (e.g., normalization, PCA,
clustering)

This registry enables compact, composable pipeline definitions (e.g.,
via character strings like `"onvps"` for "Create → Normalize →
VariableFeatures → Scale → PCA") and supports both standard
(`NormalizeData`) and SCTransform-based workflows.

## Usage

``` r
SCPreProcessStrategy
```

## Format

An object of class `environment` of length 11.

## Available Operations

- `o`:

  Create Seurat object from count matrix
  ([`SeuratObject::CreateSeuratObject`](https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html))

- `n`:

  Normalize data
  ([`Seurat::NormalizeData`](https://satijalab.org/seurat/reference/NormalizeData.html))

- `v`:

  Find variable features
  ([`Seurat::FindVariableFeatures`](https://satijalab.org/seurat/reference/FindVariableFeatures.html))

- `s`:

  Scale data
  ([`Seurat::ScaleData`](https://satijalab.org/seurat/reference/ScaleData.html))

- `p`:

  Run PCA
  ([`Seurat::RunPCA`](https://satijalab.org/seurat/reference/RunPCA.html))

- `c`:

  Find clusters
  ([`Seurat::FindClusters`](https://satijalab.org/seurat/reference/FindClusters.html))

- `e`:

  Find neighbors
  ([`Seurat::FindNeighbors`](https://satijalab.org/seurat/reference/FindNeighbors.html))

- `t`:

  Run t-SNE
  ([`Seurat::RunTSNE`](https://satijalab.org/seurat/reference/RunTSNE.html))

- `u`:

  Run UMAP
  ([`Seurat::RunUMAP`](https://satijalab.org/seurat/reference/RunUMAP.html))

- `r`:

  SCTransform
  ([`Seurat::SCTransform`](https://satijalab.org/seurat/reference/SCTransform.html))

- `i`:

  Integrate layers
  ([`Seurat::IntegrateLayers`](https://satijalab.org/seurat/reference/IntegrateLayers.html))

## Usage Notes

- These functions are **not intended for direct interactive use**. They
  are internal building blocks for workflow engines (i.e.,
  `SCPreProcess`).

- You can access any operation via `SCPreProcessStrategy$letter`, but
  doing so bypasses pipeline validation and error handling.

- To add more operations, use
  [`RegisterSeuratMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterSeuratMethod.md)

## See also

Other single_cell_preprocess:
[`FilterTumorCell()`](https://wanglabcsu.github.io/sigbridger/reference/FilterTumorCell.md),
[`FindRobustElbow()`](https://wanglabcsu.github.io/sigbridger/reference/FindRobustElbow.md),
[`Pattern2Colname()`](https://wanglabcsu.github.io/sigbridger/reference/Pattern2Colname.md),
[`QCPatternDetect()`](https://wanglabcsu.github.io/sigbridger/reference/QCPatternDetect.md),
[`RegisterSeuratMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterSeuratMethod.md),
[`SCAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/SCAnnotate.md),
[`SCIntegrate()`](https://wanglabcsu.github.io/sigbridger/reference/SCIntegrate.md),
[`SCPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcess.md),
[`compatible_with_3.0.2()`](https://wanglabcsu.github.io/sigbridger/reference/compatible_with_3.0.2.md)
