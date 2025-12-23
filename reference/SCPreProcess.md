# Single-Cell RNA-seq Preprocessing Pipeline

A generic function for standardized preprocessing of single-cell RNA-seq
data from multiple sources. Handles data.frame/matrix, AnnData, and
Seurat inputs with tumor cell filtering. Implements a complete analysis
pipeline from raw data to clustered embeddings.

## Usage

``` r
SCPreProcess(sc, ...)

# Default S3 method
SCPreProcess(
  sc,
  meta_data = NULL,
  column2only_tumor = NULL,
  project = "SC_Screening_Proj",
  min_cells = 400L,
  min_features = 0L,
  quality_control = TRUE,
  quality_control.pattern = c("^MT-"),
  data_filter = TRUE,
  data_filter.thresh = list(nFeature_RNA_thresh = c(200L, 6000L), percent.mt = 20L,
    percent.rp = 60L),
  normalization_method = "LogNormalize",
  scale_factor = 10000L,
  scale_features = NULL,
  selection_method = "vst",
  resolution = 0.6,
  dims = NULL,
  ...
)

# S3 method for class 'R6'
SCPreProcess(
  sc,
  meta_data = NULL,
  column2only_tumor = NULL,
  project = "SC_Screening_Proj",
  min_cells = 400L,
  min_features = 0L,
  quality_control = TRUE,
  quality_control.pattern = c("^MT-"),
  data_filter = TRUE,
  data_filter.thresh = list(nFeature_RNA_thresh = c(200L, 6000L), percent.mt = 20L,
    percent.rp = 60L),
  normalization_method = "LogNormalize",
  scale_factor = 10000L,
  scale_features = NULL,
  selection_method = "vst",
  resolution = 0.6,
  dims = NULL,
  ...
)

# S3 method for class 'Seurat'
SCPreProcess(sc, column2only_tumor = NULL, ...)
```

## Arguments

- sc:

  Input data, one of:

  - `data.frame/matrix/dgCMatrix`: Raw count matrix (features x cells)

  - `R6`: Python AnnData object, obtained via package `anndata` or
    `anndataR`

  - `Seurat`: Preprocessed Seurat object

- ...:

  Additional arguments passed to specific methods. Currently supports:

  - `verbose`: Logical indicating whether to print progress messages.
    Defaults to `TRUE`.

  - `dims_Neighbors`: Dimensions to use for `FindNeighbors`. Defaults to
    `NULL`, using `dims`.

  - `dims_TSNE`: Dimensions to use for `RunTSNE`. Defaults to `NULL`,
    using `dims`.

  - `dims_UMAP`: Dimensions to use for `RunUMAP`. Defaults to `NULL`,
    using `dims`.

- meta_data:

  A data.frame containing metadata for each cell. It will be added to
  the Seurat object as `@meta.data`. If sc is an anndata object, `obs`
  will be automatically used.

- column2only_tumor:

  A character of column names in `meta_data`, used to filter the Seurat
  object to only tumor cells. If `NULL`, no filtering is performed.

- project:

  A character of project name, used to name the Seurat object.

- min_cells:

  Minimum number of cells that must express a feature for it to be
  included in the analysis. Defaults to `400L`.

- min_features:

  Minimum number of features that must be detected in a cell for it to
  be included in the analysis. Defaults to `0L`.

- quality_control:

  Logical, for identification of the proportion of mitochondrial genes,
  ribosomal protein genes, or other types of genes (without filtering),
  the results will be stored in `meta.data`. Defaults to `TRUE`.

- quality_control.pattern:

  A vector or character containing regex pattern(s) to identify
  mitochondrial genes, ribosomal protein genes, or other unwanted genes,
  as well as combinations of these genes. Customized patterns are
  supported. Defaults to `"^MT-"`.

- data_filter:

  Logical indicating whether to filter cells based on quality metrics.
  Defaults to `TRUE`.

- data_filter.thresh:

  A list containing filtering thresholds for different quality metrics:

  - `nFeature_RNA_thresh`: Numeric vector of length 2 specifying the
    minimum and maximum number of features per cell. Defaults to
    `c(200L, 6000L)`

  - `percent.mt`: Maximum mitochondrial percentage allowed. Defaults to
    `20L`

  - `percent.rp`: Maximum ribosomal protein percentage allowed. Not used
    in default unless ribosomal protein genes filter pattern (like
    `^RP[LS]`) is added to `quality_control.pattern`

- normalization_method:

  Method for normalization: "LogNormalize", "CLR", or "RC". Defaults to
  `"LogNormalize"`.

- scale_factor:

  Scaling factor for normalization. Defaults to `10000L`.

- scale_features:

  Features to use for scaling. If NULL, uses all variable features.
  Defaults to `NULL`.

- selection_method:

  Method for variable feature selection: "vst", "mvp", or "disp".
  Defaults to `"vst"`.

- resolution:

  Resolution parameter for clustering. Higher values lead to more
  clusters. Defaults to `0.6`.

- dims:

  Dimensions to use for clustering and dimensionality reduction. If
  NULL, automatically determined by elbow method. Defaults to `NULL`.

## Value

A Seurat object containing:

- Normalized and scaled expression data

- Variable features identified by selection method

- PCA, t-SNE, and UMAP dimensionality reductions

- Cluster identities at specified resolution

- Quality control metrics in `@meta.data`

- When tumor cells filtered: original dimensions in `@misc$raw_dim`

- Final dimensions in `@misc$self_dim`

- Quality control column names in `@misc$qc_colnames`

## Details

**Quality Control Patterns:** The function supports flexible pattern
matching for quality control, for example:

- `"^MT-"` - Mitochondrial genes (default)

- `"^RP\[LS\]"` - Ribosomal protein genes

- `"^\[rt\]rna"` - rRNA and tRNA genes

- Custom patterns using regular expressions

- Combined patterns: `"^MT-|^RP\[LS\]"` for both mitochondrial and
  ribosomal genes

**Flexible Filtering:** The filtering system dynamically adapts to
detected quality control patterns:

- Column names are automatically generated from patterns

- Multiple thresholds can be specified in `data_filter.thresh`

- Use `SigBridgeR:::Pattern2Colname` to determine correct column names
  for custom patterns if still confused

## See also

[`CreateSeuratObject`](https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html),
[`NormalizeData`](https://satijalab.org/seurat/reference/NormalizeData.html),
[`ScaleData`](https://satijalab.org/seurat/reference/ScaleData.html),
[`FindVariableFeatures`](https://satijalab.org/seurat/reference/FindVariableFeatures.html),
[`RunPCA`](https://satijalab.org/seurat/reference/RunPCA.html),
[`RunTSNE`](https://satijalab.org/seurat/reference/RunTSNE.html),
[`RunUMAP`](https://satijalab.org/seurat/reference/RunUMAP.html),
[`FindNeighbors`](https://satijalab.org/seurat/reference/FindNeighbors.html),
[`FindClusters`](https://satijalab.org/seurat/reference/FindClusters.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# Example with matrix input
counts_matrix <- matrix(rpois(1000, 5), nrow = 100, ncol = 10)
rownames(counts_matrix) <- paste0("Gene", 1:100)
colnames(counts_matrix) <- paste0("Cell", 1:10)

seurat_obj <- SCPreProcess(
  sc = counts_matrix,
  project = "TestProject",
  min_features = 50,
  resolution = 0.8
)

# Example with tumor cell filtering
metadata <- data.frame(
  cell_type = c(rep("Tumor", 5), rep("Normal", 5)),
  row.names = paste0("Cell", 1:10)
)

tumor_seurat <- SCPreProcess(
  sc = counts_matrix,
  meta_data = metadata,
  column2only_tumor = "cell_type",
  project = "TumorAnalysis"
)
} # }
```
