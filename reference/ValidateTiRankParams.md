# Validate and Prepare DoTiRank Parameters

Internal helper that handles package installation checks, mode
resolution, and default-value resolution for
[`DoTiRank()`](https://wanglabcsu.github.io/sigbridger/reference/DoTiRank.md).

## Usage

``` r
ValidateTiRankParams(
  label_type,
  phenotype_class,
  tirank_params,
  save_path,
  load_cache,
  ...
)
```

## Arguments

- label_type, phenotype_class, tirank_params:

  Forwarded from
  [`DoTiRank()`](https://wanglabcsu.github.io/sigbridger/reference/DoTiRank.md).

- save_path, load_cache:

  Forwarded from
  [`DoTiRank()`](https://wanglabcsu.github.io/sigbridger/reference/DoTiRank.md).

- ...:

  Additional dots forwarded from
  [`DoTiRank()`](https://wanglabcsu.github.io/sigbridger/reference/DoTiRank.md).

## Value

A named list with elements: `mode`, `verbose`, `seed`, `assay`,
`load_cache`, `save_cache`, `cache_config`, `dots`.

## See also

Other TiRank:
[`DoTiRank()`](https://wanglabcsu.github.io/sigbridger/reference/DoTiRank.md)
