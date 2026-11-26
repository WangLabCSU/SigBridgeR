# Validate and Prepare DoSCIPAC Parameters

Internal helper that handles package installation checks, input
validation, family resolution, and default-value resolution for
[`DoSCIPAC()`](https://wanglabcsu.github.io/sigbridger/reference/DoSCIPAC.md).

## Usage

``` r
ValidateSCIPACParams(
  label_type,
  phenotype_class,
  hvg,
  do_pca_sc,
  n_pc,
  sc_batch_col,
  resolution,
  ela_net_alpha,
  bt_size,
  ncore,
  ci_alpha,
  nfold,
  sc_data,
  ...
)
```

## Arguments

- label_type, phenotype_class:

  Forwarded from
  [`DoSCIPAC()`](https://wanglabcsu.github.io/sigbridger/reference/DoSCIPAC.md).

- hvg, do_pca_sc, n_pc, sc_batch_col, resolution, ela_net_alpha,
  bt_size:

  Forwarded from
  [`DoSCIPAC()`](https://wanglabcsu.github.io/sigbridger/reference/DoSCIPAC.md).

- ncore, ci_alpha, nfold, sc_data:

  Forwarded from
  [`DoSCIPAC()`](https://wanglabcsu.github.io/sigbridger/reference/DoSCIPAC.md).

- ...:

  Additional dots forwarded from
  [`DoSCIPAC()`](https://wanglabcsu.github.io/sigbridger/reference/DoSCIPAC.md).

## Value

A named list with elements: `phenotype_class`, `family`, `verbose`,
`seed`, `assay`, `cache_config`.
