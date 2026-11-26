# Validate and Prepare DoSIDISH Parameters

Internal helper that handles package installation checks, input
validation, deprecated argument warnings, and default-value resolution
for
[`DoSIDISH()`](https://wanglabcsu.github.io/sigbridger/reference/DoSIDISH.md).

## Usage

``` r
ValidateSIDISHParams(
  matched_bulk,
  sc_data,
  phenotype,
  label_type,
  phenotype_class = "survival",
  sidish_params,
  ...
)
```

## Arguments

- matched_bulk, sc_data, phenotype, label_type, phenotype_class,
  sidish_params:

  Forwarded from
  [`DoSIDISH()`](https://wanglabcsu.github.io/sigbridger/reference/DoSIDISH.md).

- ...:

  Additional dots forwarded from
  [`DoSIDISH()`](https://wanglabcsu.github.io/sigbridger/reference/DoSIDISH.md).

## Value

A named list with elements: `phenotype_class`, `verbose`, `seed`,
`assay`, `env_params`, `sidish_params`.

## See also

Other SIDISH:
[`DoSIDISH()`](https://wanglabcsu.github.io/sigbridger/reference/DoSIDISH.md),
[`SIDISHEnvSet()`](https://wanglabcsu.github.io/sigbridger/reference/SIDISHEnvSet.md)
