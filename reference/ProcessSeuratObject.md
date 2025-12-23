# Process a Seurat object (internal)

Normalize, find variable features, scale, and run PCA

## Usage

``` r
ProcessSeuratObject(
  obj,
  normalization_method = "LogNormalize",
  scale_factor = 10000,
  scale_features = NULL,
  selection_method = "vst",
  verbose = TRUE
)
```

## Arguments

- obj:

  Seurat object

- normalization_method:

  Normalization method ("LogNormalize", "CLR", or "RC")

- scale_factor:

  Scaling factor for normalization

- scale_features:

  Features to scale

- selection_method:

  Variable feature selection method ("vst", "mvp", or "disp")

- verbose:

  Print progress messages

## Value

Seurat object

## See also

Other single_cell_preprocess:
[`FilterTumorCell()`](https://wanglabcsu.github.io/sigbridger/reference/FilterTumorCell.md),
[`FindRobustElbow()`](https://wanglabcsu.github.io/sigbridger/reference/FindRobustElbow.md),
[`Pattern2Colname()`](https://wanglabcsu.github.io/sigbridger/reference/Pattern2Colname.md),
[`QCPatternDetect()`](https://wanglabcsu.github.io/sigbridger/reference/QCPatternDetect.md)
