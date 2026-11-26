# Validate and Prepare DoscPAS Parameters

Internal helper that handles package installation checks, deprecated
argument warnings, input validation, and default-value resolution for
[`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md).

## Usage

``` r
ValidatescPASParams(
  phenotype,
  label_type,
  phenotype_class,
  assay,
  imputation,
  imputation_method = c("KNN", "ALRA"),
  nfeature,
  alpha,
  cutoff,
  network_class = c("SC", "bulk"),
  family = lifecycle::deprecated(),
  permutation_times,
  FDR_threshold,
  independent,
  ...
)
```

## Arguments

- phenotype, label_type, phenotype_class, assay:

  Forwarded from
  [`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md).

- imputation, imputation_method, nfeature, alpha, cutoff, network_class:

  Forwarded from
  [`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md).

- family:

  Deprecated argument forwarded from
  [`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md).

- permutation_times, FDR_threshold, independent:

  Forwarded from
  [`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md).

- ...:

  Additional dots forwarded from
  [`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md).

## Value

A named list with elements: `label_type`, `family`, `imputation_method`,
`network_class`, `verbose`, `seed`.

## See also

Other scPAS:
[`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md)
