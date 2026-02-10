# Calculate Percentage of Features Matching Patterns

This function calculates the percentage of counts coming from features
matching specified patterns (e.g., mitochondrial genes, ribosomal genes)
and adds them as metadata columns to the Seurat object.

## Usage

``` r
QCPatternDetect(
  obj,
  pattern = c("^MT-", "^mt-", "^RP[SL]", "^MT-|^RP[SL]"),
  verbose = TRUE,
  ...
)
```

## Arguments

- obj:

  A seurat object.

- pattern:

  A character vector or list containing regex patterns to identify
  mitochondrial genes, ribosomal protein genes, or other unwanted genes,
  as well as combinations of these genes. Customized patterns are
  supported.

- verbose:

  logical, whether to print progress messages

- ...:

  Additional arguments passed to
  [`PercentageFeatureSet`](https://satijalab.org/seurat/reference/PercentageFeatureSet.html)

## Details

The function automatically generates friendly column names based on the
patterns:

- "mt" for mitochondrial patterns

- "rp" for ribosomal patterns

- "rrna" for ribosomal RNA patterns

- For combined patterns (using \|), creates names like "mt_rp"

- For other patterns, creates cleaned lowercase names

## See also

Other single_cell_preprocess:
[`FilterTumorCell()`](https://wanglabcsu.github.io/sigbridger/reference/FilterTumorCell.md),
[`FindRobustElbow()`](https://wanglabcsu.github.io/sigbridger/reference/FindRobustElbow.md),
[`Pattern2Colname()`](https://wanglabcsu.github.io/sigbridger/reference/Pattern2Colname.md),
[`RegisterSeuratMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterSeuratMethod.md),
[`SCAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/SCAnnotate.md),
[`SCIntegrate()`](https://wanglabcsu.github.io/sigbridger/reference/SCIntegrate.md),
[`SCPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcess.md),
[`SCPreProcessStrategy`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcessStrategy.md),
[`compatible_with_3.0.2()`](https://wanglabcsu.github.io/sigbridger/reference/compatible_with_3.0.2.md)
