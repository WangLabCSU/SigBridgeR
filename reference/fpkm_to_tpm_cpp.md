# Convert FPKM to TPM from a beachmat Matrix

Converts every column to TPM and reports validation problems through
`status` instead of raising an error.

## Usage

``` r
fpkm_to_tpm_cpp(initialized_fpkm, na_as_zero = TRUE)
```

## Arguments

- initialized_fpkm:

  A beachmat-compatible FPKM matrix, genes in rows.

- na_as_zero:

  `logical`. Whether NA values are treated as zero. Default: `TRUE`.

## Value

A list of six elements: `tpm` (same-dimension numeric matrix on success,
`NULL` otherwise), `status` (integer: 0 success, 1 NA not allowed, 2
negative value, 3 Inf, 4 column-sum overflow, 5 scaling overflow),
`row`/`col` (1-based first problem position), `na_count` and
`zero_cols`.
