# Integrate Single-Cell Datasets

Unified interface for integrating heterogeneous single-cell datasets.
Supports `Seurat` objects and matrices (via gene-wise aggregation and
column concatenation).

## Usage

``` r
SCIntegrate(
  ...,
  pipeline = "nsvpiectu",
  add.cell.ids = NULL,
  collapse = FALSE,
  merge.data = TRUE,
  merge.dr = FALSE,
  project = "Integrated Seurat"
)

SCIntegrate.matrix(..., .quos = NULL)

SCIntegrate.Matrix(..., .quos = NULL)

SCIntegrate.Seurat(
  ...,
  pipeline = "nsvpiectu",
  add.cell.ids = NULL,
  collapse = FALSE,
  merge.data = TRUE,
  merge.dr = FALSE,
  project = "Integrated Seurat",
  .quos = NULL
)
```

## Arguments

- ...:

  For matrix method: matrices to integrate, optionally named. For Seurat
  method: Seurat objects to integrate, followed by additional parameters
  passed to integration functions, Seurat objects and parameters are
  automatically filtered and passed to appropriate functions.

- pipeline:

  Character string specifying processing steps for Seurat integration.
  Each letter represents a step: `"n"` = NormalizeData, `"s"` =
  ScaleData, `"v"` = FindVariableFeatures, `"p"` = RunPCA, `"i"` =
  IntegrateLayers, `"e"` = FindNeighbors (because `n` was used), `"c"` =
  FindClusters, `"t"` = RunTSNE, `"u"` = RunUMAP, `"r"` = SCTransform
  (because `t` was used),. Default: `"nsvpiectu"`.

- add.cell.ids:

  Character vector of prefixes for cell IDs when merging Seurat objects.
  Auto-inferred from argument names if not provided. (See Examples for
  details)

- collapse:

  If TRUE, merge layers of the same name together; if FALSE, appends
  labels to the layer name

- merge.data:

  Merge the data slots instead of just merging the counts (which
  requires renormalization); this is recommended if the same
  normalization approach was applied to all objects

- merge.dr:

  Choose how to handle merging dimensional reductions: - “TRUE”: merge
  dimensional reductions with the same name across objects; dimensional
  reductions with different names are added as-is - “NA”: keep
  dimensional reductions from separate objects separate; will append the
  project name for duplicate reduction names - “FALSE”: do not add
  dimensional reductions

- project:

  Project name for the Seurat object

- .quos:

  Ignore it. Please do not pass any value

## Value

- **Matrix method**: A matrix with all genes (union) and concatenated
  cells. Missing values filled with `NA`.

- **Seurat method**: Integrated `Seurat` object with processed layers
  according to `pipeline`.

## Matrix Method Details

- Duplicate genes (e.g., due to symbol aliasing) are resolved via
  [`AggregateDups`](https://wanglabcsu.github.io/sigbridger/reference/aggregate-dups.md)
  (default: sum).

- Dataset identifiers are auto-inferred from argument names or object
  names (e.g., `SCIntegrate(matA, B = matB)` → prefixes `"matA"`,
  `"B"`).

## Examples

``` r
if (FALSE) { # \dontrun{
# Matrix integration
mat1 <- matrix(rpois(100, 5), nrow = 20,
               dimnames = list(paste0("G", 1:20), paste0("C", 1:5)))
mat2 <- matrix(rpois(120, 6), nrow = 20,
               dimnames = list(paste0("G", 11:30), paste0("C", 1:6)))
integrated <- SCIntegrate(mat1, mat2)
dim(integrated)  # 20 genes × 11 cells

# Add customed prefixes to colnames
integrated <- SCIntegrate(A = mat1, B = mat2)
# > colnames(integrated)
# [1] "A_C1" "A_C2" "A_C3" "A_C4" "A_C5" "B_C1" "B_C2" "B_C3" "B_C4" "B_C5" "B_C6"

# Seurat integration with custom pipeline
seu1 <- CreateSeuratObject(mat1)
seu2 <- CreateSeuratObject(mat2)
integrated_seu <- SCIntegrate(seu1, seu2, pipeline = "nsfpi")
} # }

if (FALSE) { # \dontrun{
mat1 <- matrix(rpois(100, 5), nrow = 20, dimnames = list(paste0("G", 1:20), paste0("C", 1:5)))
mat2 <- matrix(rpois(120, 6), nrow = 20, dimnames = list(paste0("G", 1:20), paste0("C", 1:6)))
integrated_mat <- SCIntegrate(mat1, mat2)
dim(integrated_mat)  # 20 genes × 11 cells
head(colnames(integrated_mat))
} # }
```
