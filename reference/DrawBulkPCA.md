# Draw PCA Plot for Bulk Expression Data

Performs PCA on transposed bulk expression data and visualizes the first
two principal components as a scatter plot with marginal density plots.
Optionally overlays batch information as point shapes (up to 4 batches).

The plot is composed via
[patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html)
from:

- A top marginal density plot for PC1.

- A right marginal density plot for PC2.

- A central scatter plot with group-colored ellipses (via
  [`ggforce::geom_mark_ellipse()`](https://ggforce.data-imaginist.com/reference/geom_mark_ellipse.html)).

## Usage

``` r
DrawBulkPCA(bulk, group, batch = NULL, show_plot = TRUE, ...)
```

## Arguments

- bulk:

  A numeric matrix or data frame of bulk expression data (genes x
  samples). Column names are used as sample identifiers.

- group:

  A vector of group labels, one per sample (column of `bulk`). Used for
  coloring points and ellipses.

- batch:

  An optional vector of batch labels, one per sample. When provided,
  points are additionally shaped by batch (up to 4 batches supported).
  Default: `NULL`.

- show_plot:

  Logical. If `TRUE` (default), the plot is printed immediately.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the combined `patchwork`/`ggplot` object.

## See also

Other input_preprocess:
[`PhenoMap()`](https://wanglabcsu.github.io/sigbridger/reference/PhenoMap.md),
[`PhenoPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/PhenoPreProcess.md),
[`SCPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcess.md),
[`SymbolConvert()`](https://wanglabcsu.github.io/sigbridger/reference/SymbolConvert.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Small example matrix: 10 genes x 8 samples
set.seed(123)
bulk <- matrix(
  rnorm(100 * 80, mean = 10, sd = 2),
  nrow = 100, ncol = 80,
  dimnames = list(
    paste0("Gene", 1:100),
    paste0("Sample", 1:80)
  )
)

group <- rep(c("Control", "Treatment"), each = 40)

# Basic PCA plot without batch
DrawBulkPCA(bulk = bulk, group = group)

# PCA plot with batch shapes
batch <- rep(c("Batch1", "Batch2"), times = 40)
DrawBulkPCA(bulk = bulk, group = group, batch = batch)
} # }
```
