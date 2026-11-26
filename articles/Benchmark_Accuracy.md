# Benchmarking Screening Methods: Classification Accuracy

``` r

library(SigBridgeR)
library(Seurat)
```

## Overview

[`BenchAcc()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAcc.md)
compute the classification accuracy of screening methods according to a
known ground truth

It treats `"Positive"` as the positive class and every other supported
label (`"Negative"`, `"Neutral"`, `"Other"`) as the negative class. For
each screening method found in the object, it computes a full set of
binary classification metrics, so the methods can be ranked or compared
on exactly the same cells.

The function accepts three kinds of input:

- a `Seurat` object whose `meta.data` contains one column per screening
  method;
- a `ScreenMethodResult`, which is unwrapped to its `scRNA_data` slot;
- a `BenchAccConfig`, which runs every configured screening method first
  and then evaluates the merged labels in one pass.

## A small example

The example below uses a tiny object with six cells, a ground-truth
column and two simulated screening results. In a real analysis the
`Scissor` and `scPAS` columns would be produced by
[`Screen()`](https://wanglabcsu.github.io/sigbridger/reference/Screen.md)
or
[`MergeResult()`](https://wanglabcsu.github.io/sigbridger/reference/MergeResult.md).

``` r

counts <- Matrix::Matrix(
  matrix(
    seq_len(24L),
    nrow = 4L,
    dimnames = list(paste0("Gene", 1:4), paste0("Cell", 1:6))
  ),
  sparse = TRUE
)

sc <- CreateSeuratObject(counts = counts)

sc$truth <- c(
  "Positive",
  "Negative",
  "Positive",
  "Other",
  "Positive",
  "Negative"
)
sc$Scissor <- c(
  "Positive",
  "Negative",
  "Positive",
  "Other",
  "Negative",
  "Negative"
)
sc$scPAS <- c(
  "Positive",
  "Positive",
  "Negative",
  "Other",
  "Negative",
  "Negative"
)

metrics <- BenchAcc(sc, bench_config = "truth", verbose = FALSE)

round(do.call(rbind, metrics), 3)
#>           TPR   FPR   TNR   FNR  F1 Accuracy Precision Recall
#> Scissor 0.667 0.000 1.000 0.333 0.8    0.833       1.0  0.667
#> scPAS   0.333 0.333 0.667 0.667 0.4    0.500       0.5  0.333
```

Each element of the returned list is a one-row matrix whose row name is
the screening method. The same numbers are easier to compare after
binding them:

``` r

metrics_dt <- data.table::rbindlist(
  lapply(metrics, as.data.frame),
  idcol = "method"
)

metrics_dt[, .(method, TPR, FPR, F1, Accuracy, Precision, Recall)]
#>     method       TPR       FPR    F1  Accuracy Precision    Recall
#>     <char>     <num>     <num> <num>     <num>     <num>     <num>
#> 1: Scissor 0.6666667 0.0000000   0.8 0.8333333       1.0 0.6666667
#> 2:   scPAS 0.3333333 0.3333333   0.4 0.5000000       0.5 0.3333333
```

## Specifying the ground truth

`bench_config` can be given in two equivalent forms:

``` r

# * a length-one column name in sc@meta.data
from_column <- BenchAcc(sc, bench_config = "truth", verbose = FALSE)

# * a per-cell label vector, one label per cell
from_vector <- BenchAcc(sc, bench_config = sc$truth, verbose = FALSE)

identical(from_column, from_vector)
#> [1] TRUE
```

The labels must contain `"Positive"`, and the vector form must have one
element per cell. The function stops with an informative error when the
column is missing, when the labels have no positive class, or when the
vector length does not match the number of cells:

``` r

BenchAcc(sc, bench_config = "missing_column", verbose = FALSE)
#> Error in `BenchAcc()`:
#> ✖ "missing_column" is not a valid column name in obj@meta.data
```

## Reading the metrics

With $`P`$ positive and $`N`$ negative ground-truth cells, and $`TP`$,
$`FP`$, $`TN`$, $`FN`$ the usual confusion-matrix counts,
[`BenchAcc()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAcc.md)
reports:

| Metric | Definition | Meaning |
|----|----|----|
| `TPR` / `Recall` | $`TP / (TP + FN)`$ | Fraction of truly positive cells recovered |
| `FPR` | $`FP / (FP + TN)`$ | Fraction of negative cells incorrectly called positive |
| `TNR` | $`TN / (TN + FP)`$ | Fraction of negative cells correctly left out |
| `FNR` | $`FN / (FN + TP)`$ | Fraction of positive cells missed |
| `Precision` | $`TP / (TP + FP)`$ | Fraction of called cells that are truly positive |
| `F1` | $`2TP / (2TP + FP + FN)`$ | Harmonic mean of precision and recall |
| `Accuracy` | $`(TP + TN) / (P + N)`$ | Overall agreement with the ground truth |

Because screening methods usually return only a small positive set,
`Accuracy` tends to be high and similar across methods. `TPR`/`Recall`,
`Precision` and `F1` are usually more informative when comparing
methods. A method that selects very few cells can have perfect precision
but low recall, while a permissive method can have high recall but many
false positives.

[`BenchAcc()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAcc.md)
only scores screening method columns. Other metadata columns, including
the ground-truth column itself, are ignored. Method names are matched
case-insensitively against the methods registered in `ScreenStrategy`:

``` r

sc$SCISSOR <- sc$Scissor
sc$Scissor <- NULL

names(BenchAcc(sc, bench_config = "truth", verbose = FALSE))
#> [1] "scPAS"   "SCISSOR"
```

## Benchmarking the screening step

If the screening has not been run yet, wrap the inputs in a
`BenchAccConfig`. The configuration stores the bulk data, the
single-cell object, the phenotype and the methods to run.
[`BenchAcc()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAcc.md)
then screens the data with each method and evaluates the accumulated
labels:

``` r

config <- BenchAccConfig(
  matched_bulk = bulk,
  sc_data = sc,
  phenotype = phenotype,
  label_type = "{method_name}",
  phenotype_class = "binary",
  screen_method = c("Scissor", "scPAS"),
  method_param = MethodParamSet(
    Scissor_param = list(alpha = 0.05),
    scPAS_param = list(alpha = 0.01)
  ),
  bench_config = "truth"
)

acc <- BenchAcc(config, verbose = FALSE)
```

The `bench_config` of a `BenchAccConfig` is validated when the object is
created, so an invalid ground truth fails early instead of after a long
screening run. `method_param` is the place for method-specific
arguments; see
[`?MethodParamSet`](https://wanglabcsu.github.io/sigbridger/reference/MethodParamSet.md)
for the full list of fields.

## Practical notes

- **Use the same cells for every method.**
  [`BenchAcc()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAcc.md)
  compares the labels cell by cell. Merge the screening results first
  (for example with
  [`MergeResult()`](https://wanglabcsu.github.io/sigbridger/reference/MergeResult.md))
  so that every method shares the same cell names and order.
- **Use an independent ground truth.** A label derived from the same
  bulk phenotype that drove the screening is not an independent
  validation; use marker-based annotation, an external cell-type
  reference or held-out experimental labels when possible.
- **Report more than accuracy.** `TPR`, `Precision` and `F1` describe
  the trade-off between missing positive cells and including negative
  ones, which is usually the relevant question for phenotype-associated
  cell screening.
- **Only the four supported labels are accepted.** Predicted labels must
  be one of `"Positive"`, `"Negative"`, `"Neutral"` or `"Other"`; typos
  and `NA` abort with the offending column and row.

## See also

- [`BenchGSEA()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEA.md)
  for comparing the functional signal of the selected cells.
- [`BenchGeneScore()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScore.md)
  for comparing pathway activity at the single-cell level.
- [`BenchMemory()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemory.md)
  for comparing run time and memory.
- [`BenchAccConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAccConfig.md),
  [`MethodParamSet()`](https://wanglabcsu.github.io/sigbridger/reference/MethodParamSet.md)
  and
  [`?Screen`](https://wanglabcsu.github.io/sigbridger/reference/Screen.md)
  for the configuration and screening interfaces.
