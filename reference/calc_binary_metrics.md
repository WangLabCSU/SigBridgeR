# Calculate Binary Classification Metrics for Screening Labels

Treats `"Positive"` as the positive class and computes accuracy metrics
for each column of predicted labels.

## Usage

``` r
calc_binary_metrics(metrics, bench_label, n_threads = 1L)
```

## Arguments

- metrics:

  A data frame of predicted-label columns, one per screening method.

- bench_label:

  Character vector of ground-truth labels, one per row of `metrics`.

- n_threads:

  `integer`. Number of OpenMP threads. Default: `1L`.

## Value

A named list with one element per column of `metrics`. Each element is a
1-row numeric matrix with columns `TPR`, `FPR`, `TNR`, `FNR`, `F1` and
`Accuracy` and the method name as row name.
