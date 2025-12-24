# Filter Seurat object cells by QC metrics

Filters cells based on nFeature_RNA and optional QC metrics (e.g.
percent.mt, percent.rp), defined in `seurat_obj@misc$qc_colnames` (See
[`QCPatternDetect`](https://wanglabcsu.github.io/sigbridger/reference/QCPatternDetect.md)).
Only metrics with non-constant, non-all-zero values are used.

## Usage

``` r
QCFilter(
  seurat_obj,
  data_filter.thresh = list(nFeature_RNA_thresh = c(200L, 6000L), nCount_RNA_thresh =
    c(500L, 50000L), percent.mt = 20L, percent.rp = 60L),
  verbose = TRUE,
  ...
)
```

## Arguments

- seurat_obj:

  A `Seurat` object.

- data_filter.thresh:

  A named list with thresholds. Default:
  `list( nFeature_RNA_thresh = c(200L, 6000L), nCount_RNA_thresh = c(500L, 25000L), percent.mt = 20L, percent.rp = 60L )`.

- verbose:

  Logical; whether to print progress messages via `cli`.

- ...:

  No use

## Value

A filtered `Seurat` object.
