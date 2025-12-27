# Cluster and reduce dimensionality of single-cell data

Performs clustering and dimensionality reduction on a single-cell object
using PCA components. Automatically determines optimal dimensions when
not specified.

## Usage

``` r
ClusterAndReduce(
  obj,
  dims = NULL,
  dims_Neighbors = NULL,
  dims_TSNE = NULL,
  dims_UMAP = NULL,
  resolution = 0.6,
  verbose = TRUE,
  ...
)
```

## Arguments

- obj:

  A single-cell seurat object containing PCA reductions

- dims:

  Vector of dimensions to use (default: `NULL`, auto-determines using
  `FindRobustElbow`)

- dims_Neighbors:

  Dimensions for neighbor calculation (default: NULL, uses `dims`)

- dims_TSNE:

  Dimensions for t-SNE (default: NULL, uses `dims`)

- dims_UMAP:

  Dimensions for UMAP (default: NULL, uses `dims`)

- resolution:

  Clustering resolution parameter (default: `0.6`)

- verbose:

  Whether to print progress messages (default: `TRUE`)

- ...:

  Additional arguments passed to downstream methods, currently not used

## Value

Modified single-cell object with clustering and dimensionality reduction
results

## Note

Automatically adjusts input dimensions if they exceed available PCA
dimensions
