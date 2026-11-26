# Validate and Prepare DoscPP Parameters

Internal helper that handles package installation checks, input
validation, and default-value resolution for
[`DoscPP()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPP.md).

## Usage

``` r
ValidatescPPParams(
  phenotype,
  label_type,
  phenotype_class,
  ref_group,
  Log2FC_cutoff,
  estimate_cutoff,
  probs,
  ...
)
```

## Arguments

- phenotype, label_type, phenotype_class:

  Forwarded from
  [`DoscPP()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPP.md).

- ref_group, Log2FC_cutoff, estimate_cutoff, probs:

  Forwarded from
  [`DoscPP()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPP.md).

- ...:

  Additional dots forwarded from
  [`DoscPP()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPP.md).

## Value

A named list with elements: `phenotype_class`, `verbose`, `parallel`,
`seed`, `assay`.

## See also

Other scPP:
[`DoscPP()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPP.md)
