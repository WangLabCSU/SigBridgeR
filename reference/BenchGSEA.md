# Benchmark Screening Methods: Gene Set Enrichment Analysis

Runs gene set enrichment analysis (GSEA) on the marker genes of every
screening method present in the input, so that the functional signals
recovered by different methods can be compared side by side.

## Usage

``` r
BenchGSEA(obj, ...)

# S3 method for class 'Seurat'
BenchGSEA(
  obj,
  group,
  pathway,
  p_value_threshold = NULL,
  ...,
  verbose = getFuncOption("verbose")
)

# S3 method for class '`SigBridgeR::ScreenMethodResult`'
BenchGSEA(obj, ...)

# S3 method for class '`SigBridgeR::BenchGSEAConfig`'
BenchGSEA(obj, ..., verbose = getFuncOption("verbose"))

# Default S3 method
BenchGSEA(obj, ...)
```

## Arguments

- obj:

  A `Seurat` object whose `meta.data` contains screening method results
  (columns named after screening methods, e.g. `"Scissor"`), a
  `ScreenMethodResult` object, or a
  [BenchGSEAConfig](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEAConfig.md)
  object.

- ...:

  Additional arguments, automatically split by formal argument names
  between
  [`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/FindMarkers.html)
  and [`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html)
  (including
  [`fgsea::fgseaMultilevel()`](https://rdrr.io/pkg/fgsea/man/fgseaMultilevel.html)
  /
  [`fgsea::fgseaSimple()`](https://rdrr.io/pkg/fgsea/man/fgseaSimple.html)).
  `ident.1`, `ident.2` and `group.by` are reserved by `BenchGSEA` and
  must not be supplied.

- group:

  A vector of length equal to the number of cells, with each element
  being the group label of the corresponding cell.

- pathway:

  A named list of gene sets, each element a character vector of gene
  symbols. Passed to
  [`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html).

- p_value_threshold:

  Numeric or `NULL`. When not `NULL`, markers are filtered by adjusted
  p-value (`p_val_adj < 0.05`) before ranking. Default: `NULL`.

- verbose:

  Logical. Whether to report progress. Default: `TRUE`.

## Value

A `data.table` of
[`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html) results
stacked across screening methods, with an extra `screenMethod` column
identifying the method and a `leadingEdgeStats` list column holding the
ranking statistics of the leading-edge genes. Methods without markers
are dropped with a warning.

## Details

For each screening method found in `meta.data`:

1.  Labels are binarised to `"Positive"` vs `"Other"`;

2.  [`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/FindMarkers.html)
    is run with `ident.1 = "Positive"`, `ident.2 = "Other"` and
    `group.by` set to the method column;

3.  Markers are ranked by `sign(avg_log2FC) * -log10(p_val_adj)`;

4.  [`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html) is run
    on the ranked list, and the ranking statistics of the leading-edge
    genes are attached as `leadingEdgeStats`.

For a
[BenchGSEAConfig](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEAConfig.md)
input, every configured screening method is first run via
[`Screen()`](https://wanglabcsu.github.io/sigbridger/reference/Screen.md)
on the same data before the GSEA comparison.

## See also

[`BenchGSEAConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEAConfig.md),
[`Screen()`](https://wanglabcsu.github.io/sigbridger/reference/Screen.md),
[`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html),
[`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/FindMarkers.html),
[`BenchAcc()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAcc.md),
[`BenchMemory()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemory.md)

Other Benchmark:
[`BenchAcc()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAcc.md),
[`BenchGeneScore()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScore.md),
[`BenchMemory()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemory.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# sc is a Seurat object already screened by e.g. Scissor and scPAS
pathways <- list(
  pathway_a = c("Gene1", "Gene2", "Gene3"),
  pathway_b = c("Gene4", "Gene5", "Gene6")
)
res <- BenchGSEA(sc, pathway = pathways)
} # }
```
