# Automatically determine optimal PCA dimensions using multiple robust methods

This function combines multiple statistical approaches to automatically
determine the optimal number of principal components (PCs) for
downstream single-cell analysis. It integrates variance-based
heuristics, elbow detection algorithms, and provides comprehensive
visualization for result validation.

## Usage

``` r
FindRobustElbow(
  obj,
  verbose = SigBridgeRUtils::getFuncOption("verbose") %||% TRUE,
  ndims = 50L
)
```

## Arguments

- obj:

  A Seurat object that has PCA computed (after `RunPCA`)

- verbose:

  Logical, if TRUE outputs detailed method results and creates
  visualization plot. If FALSE returns only the final dimension.

- ndims:

  Integer, maximum number of dimensions to consider (default: `50L`)

## Value

Integer, the recommended number of PCA dimensions for downstream
analysis

## See also

Other single_cell_preprocess:
[`FilterTumorCell()`](https://wanglabcsu.github.io/sigbridger/reference/FilterTumorCell.md),
[`Pattern2Colname()`](https://wanglabcsu.github.io/sigbridger/reference/Pattern2Colname.md),
[`QCPatternDetect()`](https://wanglabcsu.github.io/sigbridger/reference/QCPatternDetect.md),
[`SCPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcess.md),
[`compatible_with_3.0.2()`](https://wanglabcsu.github.io/sigbridger/reference/compatible_with_3.0.2.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# After running PCA on Seurat object
pbmc <- RunPCA(pbmc, npcs = 50)
optimal_dims <- FindRobustElbow(pbmc, verbose = TRUE)
pbmc <- FindNeighbors(pbmc, dims = 1:optimal_dims)
} # }
```
