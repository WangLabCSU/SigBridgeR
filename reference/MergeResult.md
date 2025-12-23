# Merge Multiple Screening Analysis Results

Combines results from multiple single-cell screening analyses (Scissor,
scPAS, scPP, or scAB) by merging their metadata and miscellaneous
information while preserving the original expression data. Performs an
inner join on cell barcodes to ensure only cells present in all inputs
are retained.

## Usage

``` r
MergeResult(..., verbose = SigBridgeRUtils::getFuncOption("verbose"))
```

## Arguments

- ...:

  Input objects to merge. Can be: - Seurat objects - Lists containing
  `scRNA_data` (Seurat objects) - Mixed combinations of the above - The
  first one will be used as base object for merging, priority given to
  first one when duplicate columns are found

- verbose:

  Logical, whether to print a message

## Value

A merged Seurat object containing:

- Expression data from the first input object

- Combined metadata from all input objects

- Miscellaneous information from all input objects

- Only cells present in all input objects (inner join)

## Processing Details

1.  Input Validation: Checks for valid Seurat objects or lists
    containing Seurat objects

2.  Metadata Extraction: Collects metadata from all objects

3.  Cell Intersection: Retains only cells present in all datasets

4.  Object Merging: Creates new Seurat object with combined metadata

5.  Miscellaneous: Adds miscellaneous information to the merged object

## Examples

``` r
if (FALSE) { # \dontrun{
# Merge mixed analysis types
combined <- MergeResult(scissor_output, scAB_output, scPP_output)

# Merge list-containing objects
merged_list <- MergeResult(list1, list2, seurat_obj)
} # }
```
