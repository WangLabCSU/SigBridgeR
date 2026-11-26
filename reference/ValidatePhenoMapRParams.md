# Validate and Prepare DoPhenoMapR Parameters

Internal helper that handles package checks, sample-id alignment,
phenotype normalization, input validation, and default-value resolution
for
[`DoPhenoMapR()`](https://wanglabcsu.github.io/sigbridger/reference/DoPhenoMapR.md).

## Usage

``` r
ValidatePhenoMapRParams(
  matched_bulk,
  phenotype,
  label_type,
  phenotype_class,
  verbose,
  seed,
  ...
)
```

## Arguments

- matched_bulk:

  Matrix or data frame of preprocessed bulk RNA-seq expression data
  (genes x samples). Column names must match names/IDs in phenotype.

- phenotype:

  Phenotype data, either:

  - Named vector, where names match columns of matched_bulk

  - Patient survival data frame, where row names match columns of
    matched_bulk, and column names are "time" and "status"

- label_type:

  Character specifying phenotype label type, e.g. "SBS1" or "time".

- phenotype_class:

  Type of phenotypic outcome. Must be consistent with input data:

  - "binary": Binary traits, e.g. case/control

  - "continuous": Continuous measurements

  - "survival": Survival information

- ...:

  Additional dots forwarded from
  [`DoPhenoMapR()`](https://wanglabcsu.github.io/sigbridger/reference/DoPhenoMapR.md).

## Value

A named list with the normalized phenotype data frame (`pheno_df`), the
derived parameters, and a `ScreenMethodConfig` (`cache_config`).

## See also

Other PhenoMapR:
[`DoPhenoMapR()`](https://wanglabcsu.github.io/sigbridger/reference/DoPhenoMapR.md)
