# Single-Cell RNA-seq Preprocessing Pipeline

A unified interface for standardized single-cell RNA-seq preprocessing.
It accepts raw counts (matrix/data.frame), AnnData (via `anndata` or
`anndataR`), or existing Seurat objects. The analysis flow is fully
customizable via a character string `pipeline` and a configuration list
`params`.

## Usage

``` r
SCPreProcess(sc, ...)

# Default S3 method
SCPreProcess(
  sc,
  ...,
  pipeline = "onsvpcetu",
  params = list(o = list(project = "SC_Screen_Proj", min.cells = 400L), n = list(), s =
    list(), v = list(), p = list(), e = list(), c = list(resolution = 0.6), t = list(), u
    = list()),
  quality_control = list(pattern = c("^MT-")),
  data_filter = list(nFeature_RNA_thresh = c(200L, 6000L), nCount_RNA_thresh = c(500L,
    50000L), percent.mt = 20L, percent.rp = 60L),
  column2only_tumor = NULL
)

# S3 method for class 'R6'
SCPreProcess(
  sc,
  ...,
  pipeline = "onsvpcetu",
  params = list(o = list(project = "SC_Screen_Proj", min.cells = 400L), n = list(), s =
    list(), v = list(), p = list(), e = list(), c = list(resolution = 0.6), t = list(), u
    = list()),
  quality_control = list(pattern = c("^MT-")),
  data_filter = list(nFeature_RNA_thresh = c(200L, 6000L), nCount_RNA_thresh = c(500L,
    50000L), percent.mt = 20L, percent.rp = 60L),
  column2only_tumor = NULL
)

# S3 method for class 'Seurat'
SCPreProcess(sc, column2only_tumor = NULL, ...)
```

## Arguments

- sc:

  Input data. Can be:

  - **Matrix/Data frame**: Raw count matrix (genes x cells).

  - **AnnData**: Python AnnData object (read via `anndata` or `anndataR`
    packages).

  - **Seurat**: A Seurat object (automatically validated and repaired if
    necessary).

- ...:

  Additional arguments for backward compatibility (mapped to `params`)
  or verbose control.

- pipeline:

  A character string defining the processing steps and order. Characters
  map to Seurat functions:

  - `'o'`: `CreateSeuratObject` (Must be the first step and cannot be
    deleted)

  - `'n'`: `NormalizeData`

  - `'s'`: `ScaleData`

  - `'v'`: `FindVariableFeatures`

  - `'p'`: `RunPCA`

  - `'e'`: `FindNeighbors` (Because "n" is used)

  - `'c'`: `FindClusters`

  - `'t'`: `RunTSNE`

  - `'u'`: `RunUMAP`

  - `'r'`: `SCTransform` (Alternative to n/s/v)

  Default is `"onsvpcetu"`.

- params:

  A named list of lists containing arguments for each pipeline step.
  Keys match the pipeline characters (e.g., `params$n` for
  `NormalizeData`). Default structure:

      list(
        o = list(project = "SC_Screen_Proj", min.cells = 400L), # do not pass `counts`
        n = list(),             # NormalizeData args
        s = list(),             # ScaleData args
        v = list(),             # FindVariableFeatures args
        c = list(resolution = 0.6),
        ...
      )

- quality_control:

  A list containing regex patterns for QC metric calculation. (See
  [QCPatternDetect](https://wanglabcsu.github.io/sigbridger/reference/QCPatternDetect.md))
  Default: `list(pattern = "^MT-")`. Detected metrics (e.g., percent.mt)
  are added to meta.data.

- data_filter:

  A list of thresholds for cell filtering. Default: `nFeature_RNA`
  (200-6000), `nCount_RNA` (500-50000), `percent.mt` (\<20),
  `percent.rp` (\<60). Only metrics detected via `quality_control` are
  filtered, i.e., nFeature_RNA, nCount_RNA and percent.mt.

- column2only_tumor:

  Optional character. Column name in metadata to filter for tumor cells
  (matches "Tumor", "Cancer", "Malignant", etc.). If `NULL`, no
  filtering is performed.

## Value

A processed Seurat object with reductions, clusters, and QC metrics.

## Details

**Pipeline Strategy:** The function parses the `pipeline` string and
executes corresponding Seurat functions in order. To use `SCTransform`,
simply change the pipeline string (e.g., `"orpetu"`) and provide
parameters in `params$r`.

**Quality Control & Filtering:** QC metrics are generated based on regex
patterns in `quality_control`. Cells are then filtered based on
thresholds in `data_filter`. Column names for filtering are
auto-generated (e.g., pattern "^MT-" -\> filter "percent.mt"). If
confused about the column name, use `SigBridgeR:::Pattern2Colname()`.

## See also

Other single_cell_preprocess:
[`FilterTumorCell()`](https://wanglabcsu.github.io/sigbridger/reference/FilterTumorCell.md),
[`FindRobustElbow()`](https://wanglabcsu.github.io/sigbridger/reference/FindRobustElbow.md),
[`Pattern2Colname()`](https://wanglabcsu.github.io/sigbridger/reference/Pattern2Colname.md),
[`QCPatternDetect()`](https://wanglabcsu.github.io/sigbridger/reference/QCPatternDetect.md),
[`compatible_with_3.0.2()`](https://wanglabcsu.github.io/sigbridger/reference/compatible_with_3.0.2.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# 1. Standard pipeline (LogNormalize -> Scale -> PCA -> UMAP)
obj <- SCPreProcess(
  sc = counts_matrix,
  pipeline = "onsvpcetu",
  params = list(c = list(resolution = 0.8))
)

# 2. SCTransform pipeline
obj_sct <- SCPreProcess(
  sc = counts_matrix,
  pipeline = "orpcu", # Create -> SCT -> PCA -> Clusters -> UMAP
  quality_control = list(pattern = c("^MT-", "^RP[LS]")),
  params = list(
    r = list(vars.to.regress = "percent.mt")
  )
)

# 3. Start from AnnData with tumor filtering
adata_object <- anndataR::read_h5ad("data.h5ad")
obj_ad <- SCPreProcess(
  sc = adata_object,
  column2only_tumor = "tissue"
)
} # }
```
