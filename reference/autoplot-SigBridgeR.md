# autoplot method for Seurat UMAP

A ggplot2 autoplot method for class ScreenMethodResult and SeuratObject

## Usage

``` r
# S3 method for class 'Seurat'
autoplot(
  object,
  reduction = "umap",
  group.by = NULL,
  label = FALSE,
  label.size = 4L,
  pt.size = NULL,
  cols = NULL,
  ...
)

# S3 method for class '`SigBridgeR::ScreenMethodResult`'
autoplot(object, ...)
```

## Arguments

- object:

  Seurat object

- reduction:

  Which dimensionality reduction to use. If not specified, first
  searches for umap, then tsne, then pca

- group.by:

  Name of one or more metadata columns to group (color) cells by (for
  example, orig.ident); pass 'ident' to group by identity class

- label:

  Whether to label the clusters

- label.size:

  Sets size of labels

- pt.size:

  Adjust point size for plotting

- cols:

  Vector of colors, each color corresponds to an identity class. This
  may also be a single character or numeric value corresponding to a
  palette as specified by
  [`brewer.pal.info`](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html).
  By default, ggplot2 assigns colors. We also include a number of
  palettes from the pals package. See
  [`DiscretePalette`](https://satijalab.org/seurat/reference/DiscretePalette.html)
  for details.

- ...:

  Extra parameters passed to `DimPlot`

## Value

ggplot object
