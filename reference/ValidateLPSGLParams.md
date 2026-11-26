# Validate and Prepare DoLP_SGL Parameters

Internal helper that handles package installation checks, deprecated
argument warnings, input validation, and default-value resolution for
[`DoLP_SGL()`](https://wanglabcsu.github.io/sigbridger/reference/DoLP_SGL.md).

## Usage

``` r
ValidateLPSGLParams(
  label_type,
  phenotype_class = c("binary", "continuous", "survival"),
  family = lifecycle::deprecated(),
  resolution,
  alpha,
  nfold,
  dge_analysis,
  ...
)
```

## Arguments

- label_type, phenotype_class, family:

  Forwarded from
  [`DoLP_SGL()`](https://wanglabcsu.github.io/sigbridger/reference/DoLP_SGL.md).

- resolution, alpha, nfold, dge_analysis:

  Forwarded from
  [`DoLP_SGL()`](https://wanglabcsu.github.io/sigbridger/reference/DoLP_SGL.md).

- ...:

  Additional dots forwarded from
  [`DoLP_SGL()`](https://wanglabcsu.github.io/sigbridger/reference/DoLP_SGL.md).

## Value

A named list with elements: `phenotype_class`, `family`, `verbose`,
`seed`, `assay`, `cache_config`.

## See also

Other LP_SGL:
[`DoLP_SGL()`](https://wanglabcsu.github.io/sigbridger/reference/DoLP_SGL.md)
