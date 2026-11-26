# Validate and Prepare DoscAB Parameters

Internal helper that handles package installation checks, input
validation, and default-value resolution for
[`DoscAB()`](https://wanglabcsu.github.io/sigbridger/reference/DoscAB.md).

## Usage

``` r
ValidatescABParams(
  phenotype,
  label_type,
  phenotype_class = c("binary", "survival"),
  alpha,
  alpha_2,
  k_max,
  cross_k,
  repeat_times,
  maxiter,
  tred,
  ...
)
```

## Arguments

- phenotype, label_type, phenotype_class:

  Forwarded from
  [`DoscAB()`](https://wanglabcsu.github.io/sigbridger/reference/DoscAB.md).

- alpha, alpha_2, k_max, cross_k, repeat_times, maxiter, tred:

  Forwarded from
  [`DoscAB()`](https://wanglabcsu.github.io/sigbridger/reference/DoscAB.md).

- ...:

  Additional dots forwarded from
  [`DoscAB()`](https://wanglabcsu.github.io/sigbridger/reference/DoscAB.md).

## Value

A named list with elements: `phenotype_class`, `verbose`, `seed`,
`parallel`, `assay`, `load_cache`, `save_cache`, `cache_config`, `dots`.

## See also

Other scAB:
[`DoscAB()`](https://wanglabcsu.github.io/sigbridger/reference/DoscAB.md)
