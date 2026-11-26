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

  A Seurat object that has PCA computed (after `RunPCA`) or a numeric
  matrix of cell embeddings (cells x PCs).

- verbose:

  Logical, if TRUE outputs detailed method results and creates
  visualization plot. If `FALSE` returns only the final dimension.

- ndims:

  Integer, maximum number of dimensions to consider (default: `50L`).
  Values exceeding the number of PCs stored in the `"pca"` reduction are
  capped at the available number of PCs.

## Value

Integer, the recommended number of PCA dimensions for downstream
analysis

## Details

The numerical work is delegated to the C++ routine
[`find_robust_elbow_cpp()`](https://wanglabcsu.github.io/sigbridger/reference/find_robust_elbow_cpp.md),
which receives the PCA cell embeddings as an S4-aware matrix view
created by
[`beachmat::initializeCpp()`](https://rdrr.io/pkg/beachmat/man/initializeCpp.html).
Columns are streamed through `tatami` (see `Rtatami.h`), so the
embeddings are never coerced into a dense base `matrix` and only the
requested `ndims` columns are touched. For each of those PCs the routine
computes the standard deviation, the explained variance, the percentage
of explained variance, and the cumulative percentage, which are then fed
into the six heuristics described below.

The recommended dimension is the second largest value among the
following candidates:

- cumulative variance above 90\\

- cumulative variance above 80\\

- last PC explaining more variance than the mean;

- last PC explaining more than two standard deviations of the variance;

- maximum second derivative of the variance curve;

- largest distance between the variance curve and the chord joining the
  first and the last PC;

- the fixed anchors 10, 15, 20, 25, 30 and 35.

The fixed anchors are independent of the data, so the raw recommendation
can exceed the number of PCs that were actually computed. The final
value is therefore capped at the number of available PCs (i.e. `ndims`).

When `verbose = TRUE` and the session is interactive, the per-method
candidates are plotted on top of the variance curve by the internal
helper
[`plot_robust_elbow()`](https://wanglabcsu.github.io/sigbridger/reference/plot_robust_elbow.md).

## See also

Other single_cell_preprocess:
[`FilterTumorCell()`](https://wanglabcsu.github.io/sigbridger/reference/FilterTumorCell.md),
[`Pattern2Colname()`](https://wanglabcsu.github.io/sigbridger/reference/Pattern2Colname.md),
[`QCPatternDetect()`](https://wanglabcsu.github.io/sigbridger/reference/QCPatternDetect.md),
[`RegisterSeuratMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterSeuratMethod.md),
[`SCAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/SCAnnotate.md),
[`SCIntegrate()`](https://wanglabcsu.github.io/sigbridger/reference/SCIntegrate.md),
[`SCPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcess.md),
[`SCPreProcessStrategy`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcessStrategy.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# After running PCA on Seurat object
pbmc <- RunPCA(pbmc, npcs = 50)
optimal_dims <- FindRobustElbow(pbmc, verbose = TRUE)
pbmc <- FindNeighbors(pbmc, dims = 1:optimal_dims)
} # }
```
