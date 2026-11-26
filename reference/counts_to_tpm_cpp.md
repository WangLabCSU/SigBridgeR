# Convert Counts to TPM from a beachmat Matrix

Convert Counts to TPM from a beachmat Matrix

## Usage

``` r
counts_to_tpm_cpp(initialized_counts, inv_length)
```

## Arguments

- initialized_counts:

  A beachmat-compatible counts matrix, genes in rows and samples in
  columns.

- inv_length:

  Numeric vector of `1 / gene_length`, aligned to the rows on the R
  side.

## Value

A numeric matrix with the same dimensions as `initialized_counts`,
holding TPM values.
