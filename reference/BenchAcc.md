# Benchmark Screening Methods: Classification Accuracy

Computes binary classification metrics for every screening method
present in the input, by comparing each method's per-cell labels against
ground-truth labels. `"Positive"` is treated as the positive class and
all other labels as negative, so the accuracy of different screening
methods can be compared on the same cells.

## Usage

``` r
BenchAcc(obj, ...)

# S3 method for class 'Seurat'
BenchAcc(
  obj,
  bench_config = character(),
  ...,
  verbose = getFuncOption("verbose"),
  nthreads = getFuncOption("nthreads")
)

# S3 method for class '`SigBridgeR::BenchAccConfig`'
BenchAcc(obj, ..., verbose = getFuncOption("verbose"))

# S3 method for class '`SigBridgeR::ScreenMethodResult`'
BenchAcc(obj, ...)

# Default S3 method
BenchAcc(obj, ...)
```

## Arguments

- obj:

  The object to evaluate, one of:

  - a `Seurat` object whose `meta.data` contains screening method
    results (columns named after screening methods, e.g. `"Scissor"`);

  - a
    [ScreenMethodResult](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodResult.md)
    object, whose `scRNA_data` slot is evaluated;

  - a
    [BenchAccConfig](https://wanglabcsu.github.io/sigbridger/reference/BenchAccConfig.md)
    object, whose `screen_method` entries are all run through
    [`Screen()`](https://wanglabcsu.github.io/sigbridger/reference/Screen.md)
    first and the accumulated labels are evaluated in a single pass.

- ...:

  Additional arguments passed between the S3 methods. Currently unused
  apart from forwarding

- bench_config:

  `character`. Ground-truth labels used as the reference: either a
  length-one column name in `obj@meta.data`, or a per-cell label vector
  of length `ncol(obj)`. Must contain a `"Positive"` label. For a
  [BenchAccConfig](https://wanglabcsu.github.io/sigbridger/reference/BenchAccConfig.md)
  input the labels are taken from the same-named slot.

- verbose:

  Logical. Whether to report progress information. Default: `TRUE`.

- nthreads:

  Integer. Number of threads used by the C++ metric computation.
  Default: `2L`.

## Value

A named `list` with one element per screening method found in
`meta.data`. Each element is a 1-row numeric matrix with columns `TPR`,
`FPR`, `TNR`, `FNR`, `F1`, `Accuracy` and so on.

## Details

### Dispatch

`BenchAcc()` dispatches on the class of `obj`: a
[BenchAccConfig](https://wanglabcsu.github.io/sigbridger/reference/BenchAccConfig.md)
is screened first, a
[ScreenMethodResult](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodResult.md)
is unwrapped to its `scRNA_data` slot, and `BenchAcc.Seurat()` performs
the computation. Unsupported classes abort with an informative error.

### Method detection

Columns of `obj@meta.data` whose names match a supported screening
method (case-insensitive) are taken as the predicted labels; other
columns are ignored. The function aborts when no such column exists.

### Ground truth

A length-one `bench_config` is looked up as a column of `obj@meta.data`;
otherwise it is used directly as a per-cell label vector. The labels
must include `"Positive"`, and the vector form must have one element per
cell.

### Metrics

With `P`, `N`, `TP`, `FP`, `TN` and `FN` denoting the numbers of
positive cells, negative cells, true/false positives and true/false
negatives

## See also

[`BenchAccConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAccConfig.md),
[ScreenMethodResult](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodResult.md),
[`Screen()`](https://wanglabcsu.github.io/sigbridger/reference/Screen.md),
[`BenchGSEA()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEA.md),
[`BenchMemory()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemory.md),
[`BenchGeneScore()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScore.md)

Other Benchmark:
[`BenchGSEA()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEA.md),
[`BenchGeneScore()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScore.md),
[`BenchMemory()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemory.md)

## Examples

``` r
counts <- matrix(
  rep(1:10, 6L),
  nrow = 10L,
  dimnames = list(paste0("Gene", 1:10), paste0("Cell", 1:6))
)
obj <- Seurat::CreateSeuratObject(counts = counts)
#> Warning: Data is of class matrix. Coercing to dgCMatrix.
obj$truth <- c("Positive", "Negative", "Positive", "Other", "Positive", "Negative")
obj$Scissor <- c("Positive", "Negative", "Positive", "Other", "Negative", "Negative")
BenchAcc(obj, bench_config = "truth", verbose = FALSE)
#> $`c("Positive", "Negative", "Positive", "Other", "Negative", "Negative")`
#>                                                                              TPR
#> c("Positive", "Negative", "Positive", "Other", "Negative", "Negative") 0.6666667
#>                                                                        FPR TNR
#> c("Positive", "Negative", "Positive", "Other", "Negative", "Negative")   0   1
#>                                                                              FNR
#> c("Positive", "Negative", "Positive", "Other", "Negative", "Negative") 0.3333333
#>                                                                         F1
#> c("Positive", "Negative", "Positive", "Other", "Negative", "Negative") 0.8
#>                                                                         Accuracy
#> c("Positive", "Negative", "Positive", "Other", "Negative", "Negative") 0.8333333
#>                                                                        Precision
#> c("Positive", "Negative", "Positive", "Other", "Negative", "Negative")         1
#>                                                                           Recall
#> c("Positive", "Negative", "Positive", "Other", "Negative", "Negative") 0.6666667
#> 
```
