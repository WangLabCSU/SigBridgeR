# Train DEGAS Model and Optionally Cache

Internal helper that sets up the Python environment, trains the DEGAS
model using bootstrap-aggregated multi-task learning, and optionally
persists the result to cache.

## Usage

``` r
TrainDEGASModel(
  matched_bulk,
  sc_data,
  phenotype,
  t_sc_mat,
  cm_genes,
  sc_data.pheno_colname,
  label_type,
  p
)
```

## Arguments

- matched_bulk, sc_data, phenotype, sc_data.pheno_colname, label_type:

  Forwarded from
  [`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md).

- p:

  Named list of validated parameters from
  [`ValidateDEGASParams()`](https://wanglabcsu.github.io/sigbridger/reference/ValidateDEGASParams.md).

## Value

The trained `ccModel1` object (list of DEGAS models).

## See also

Other DEGAS:
[`DEGASEnvSet()`](https://wanglabcsu.github.io/sigbridger/reference/DEGASEnvSet.md),
[`DEGASModelDetect()`](https://wanglabcsu.github.io/sigbridger/reference/DEGASModelDetect.md),
[`DEGASParamSet()`](https://wanglabcsu.github.io/sigbridger/reference/DEGASParamSet.md),
[`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md),
[`ValidateDEGASParams()`](https://wanglabcsu.github.io/sigbridger/reference/ValidateDEGASParams.md)
