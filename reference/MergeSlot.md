# Helper function to merge slot

Merges the slot of a Seurat object (`merged_obj`) with the slot
(`slot_type`) of other Seurat objects (`seurat_objects`). `common_cells`
is a vector of cell barcodes that are common to all Seurat objects.

## Usage

``` r
MergeSlot(slot_type, merged_obj, seurat_objects, common_cells)
```
