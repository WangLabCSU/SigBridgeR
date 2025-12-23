# Add Gene-Level Metadata to Seurat Object (Vectorized, ...-based)

Add multiple feature-level metadata (vectors or 2D tables) to a Seurat
object. Handles vectors (length-checked), matrices, data.frames,
data.tables, etc. Columns/variables with duplicated names are suffixed
(e.g., "type_1"). Gene alignment is auto-detected: rows or cols must
match nrow(seurat_obj).

## Usage

``` r
AddMetaFeature(seurat_obj, ..., assay = "RNA")
```

## Arguments

- seurat_obj:

  A Seurat object.

- ...:

  One or more metadata inputs:

  - Named/unnamed vectors (length = ngenes)

  - 2D objects (matrix, data.frame, data.table, etc.) where one
    dimension (rows or cols) has size = ngenes.

- assay:

  Assay name (default: `"RNA"`, fallback to first).

## Value

Modified `seurat_obj` (invisibly).
