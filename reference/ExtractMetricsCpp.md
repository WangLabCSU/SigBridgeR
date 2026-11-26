# Extract Mean-Variance and Dropout Metrics

Derives summary metrics from assay values and the corresponding counts.

## Usage

``` r
ExtractMetricsCpp(assay_data, counts, low_expressed_thresh = 0.2)
```

## Arguments

- assay_data:

  A beachmat-compatible matrix of assay values, genes in rows.

- counts:

  A beachmat-compatible counts matrix with the same dimensions as
  `assay_data`.

- low_expressed_thresh:

  `numeric`. Threshold for considering a gene lowly expressed. Default:
  `0.2`.

## Value

A list with `variance_mean_cor` and `mean_dropout_residual` (numeric
scalars) and `n_cells`, `n_genes` (integer scalars).
