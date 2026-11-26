# Validate and Prepare DoDEGAS Parameters

Internal helper that handles package installation checks, input
validation, default-value resolution, and cache configuration for
[`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md).

## Usage

``` r
ValidateDEGASParams(
  select_fraction,
  min_thresh,
  matched_bulk,
  sc_data,
  phenotype,
  sc_data.pheno_colname,
  label_type,
  phenotype_class,
  tmp_dir,
  env_params = lifecycle::deprecated(),
  degas_params,
  normality_test_method,
  ...
)
```

## Arguments

- select_fraction, min_thresh, matched_bulk, sc_data, phenotype:

  Forwarded from
  [`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md).

- sc_data.pheno_colname, label_type, phenotype_class, tmp_dir:

  Forwarded from
  [`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md).

- env_params, degas_params, normality_test_method:

  Forwarded from
  [`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md).

- ...:

  Additional dots forwarded from
  [`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md).

## Value

A named list with elements: `phenotype_class`, `normality_test_method`,
`verbose`, `assay`, `load_cache`, `save_cache`, `env_params`,
`degas_params`, `cache_config`, `dots`.

## See also

Other DEGAS:
[`DEGASEnvSet()`](https://wanglabcsu.github.io/sigbridger/reference/DEGASEnvSet.md),
[`DEGASModelDetect()`](https://wanglabcsu.github.io/sigbridger/reference/DEGASModelDetect.md),
[`DEGASParamSet()`](https://wanglabcsu.github.io/sigbridger/reference/DEGASParamSet.md),
[`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md),
[`TrainDEGASModel()`](https://wanglabcsu.github.io/sigbridger/reference/TrainDEGASModel.md)
