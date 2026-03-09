# Merge Multiple Screening Analysis Results

Combines results from multiple single-cell screening analyses (e.g.,
Scissor, scPAS, scPP, scAB) by merging their metadata, slots, and
miscellaneous information while preserving the original expression data
from the first input. Performs an inner join on cell barcodes to ensure
only cells present in all inputs are retained.

## Usage

``` r
MergeResult(
  ...,
  weights = NULL,
  ties.method = c("random", "first", "last"),
  verbose = SigBridgeRUtils::getFuncOption("verbose")
)
```

## Arguments

- ...:

  Input objects to merge. Can be:

  - Seurat objects

  - Lists containing `scRNA_data` (Seurat objects)

  - Mixed combinations of the above

  The first object serves as the base for merging. When duplicate
  metadata columns are found, priority is given to the first occurrence.

- weights:

  Named numeric vector specifying weights for voting when multiple
  screening methods produce categorical labels. Names must match
  screening method names (case-insensitive). If `NULL` (default), no
  voting is performed.

- ties.method:

  Character. Method for resolving ties in weighted voting:

  `"first"`

  :   Select the first category alphabetically (default)

  `"last"`

  :   Select the last category alphabetically

  `"random"`

  :   Randomly select among tied categories

- verbose:

  Logical. Whether to print progress messages. Default: inherits from
  `getOption("SigBridgeR.verbose")`.

## Value

A merged Seurat object containing:

- Expression data from the first input object

- Combined metadata from all input objects (inner join on cell IDs)

- Merged slots (assays, reductions, graphs, images) from all objects

- Combined miscellaneous information from all objects

- Optional `vote` column in metadata if `weight` is provided

## Processing Workflow

1.  **Input Validation**: Extracts Seurat objects from inputs (skips
    invalid objects with warning)

2.  **Metadata Extraction**: Converts metadata to data.table format for
    efficient merging

3.  **Cell Intersection**: Retains only cells present in all datasets
    (inner join)

4.  **Weighted Voting**: If `weight` provided, performs weighted voting
    across screening methods

5.  **Slot Merging**: Combines assays, reductions, graphs, and images
    from all objects

6.  **Misc Merging**: Aggregates miscellaneous information from all
    objects

## Weighted Voting

When `weight` is provided, the function:

1.  Identifies columns in metadata matching registered screening method
    names

2.  Extracts vote data (must contain "Positive", "Negative", "Neutral",
    or "Other")

3.  Applies
    [`WeightedVote`](https://wanglabcsu.github.io/sigbridger/reference/WeightedVote.md)
    with specified weights

4.  Adds `vote` column to merged metadata with final aggregated labels

## Notes

- Only cells present in *all* input objects are retained (inner join
  behavior)

- Duplicate metadata columns are resolved by keeping the first
  occurrence

- For integrating *heterogeneous* single-cell datasets (different
  samples), use
  [`SCIntegrate`](https://wanglabcsu.github.io/sigbridger/reference/SCIntegrate.md)
  instead

- Weights should be non-negative; higher values indicate greater
  influence in voting

## See also

[`WeightedVote`](https://wanglabcsu.github.io/sigbridger/reference/WeightedVote.md),
[`SCIntegrate`](https://wanglabcsu.github.io/sigbridger/reference/SCIntegrate.md),
[`AddMisc`](https://wanglabcsu.github.io/sigbridger/reference/AddMisc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Merge results from different screening methods
scissor_out <- RunScissor(...)
scAB_out <- RunscAB(...)
scPP_out <- RunscPP(...)

merged <- MergeResult(scissor_out, scAB_out, scPP_out)

# Check merged metadata
head(merged[[]])

# Example 2: Weighted voting across methods
weights <- c(
  Scissor = 3,    # Highest confidence
  scAB = 2,       # Medium confidence
  scPP = 1        # Lower confidence
)

merged <- MergeResult(
  scissor_out,
  scAB_out,
  scPP_out,
  weight = weights,
  ties.method = "first"
)

# View voting results
table(merged$vote)

# Example 3: Merge list-containing objects
result_list1 <- list(scRNA_data = seurat_obj1, other_info = "...")
result_list2 <- list(scRNA_data = seurat_obj2, other_info = "...")

merged <- MergeResult(result_list1, result_list2)

# Example 4: Check which cells were retained
original_cells <- ncol(seurat_obj1)
merged_cells <- ncol(merged)
cat("Retained", merged_cells, "of", original_cells, "cells\n")

# Example 5: Access merged slots
# All unique assays from all inputs
names(merged@assays)

# All unique reductions
names(merged@reductions)
} # }
```
