# Validate and Prepare DoPIPET Parameters

Internal helper that handles package installation checks, input
validation, and default-value resolution for
[`DoPIPET()`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md).

## Usage

``` r
ValidatePIPETParams(
  phenotype_class,
  only_pos_marker,
  group,
  discretize_method,
  cutoff,
  label_type,
  marker_finder,
  log2FC,
  p_adjust,
  show_log2FC,
  freq_counts,
  normalize,
  scale,
  nPerm,
  distance,
  ...
)
```

## Arguments

- phenotype_class, only_pos_marker, group:

  Forwarded from
  [`DoPIPET()`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md).

- discretize_method, cutoff, label_type, marker_finder, log2FC,
  p_adjust:

  Forwarded from
  [`DoPIPET()`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md).

- show_log2FC, freq_counts, normalize, scale, nPerm, distance:

  Forwarded from
  [`DoPIPET()`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md).

- ...:

  Additional dots forwarded from
  [`DoPIPET()`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md).

## Value

A named list with elements: `phenotype_class`, `marker_finder`,
`distance`, `verbose`, `seed`, `parallel`, `assay`, `load_cache`,
`save_cache`, `cache_config`, `dots`.

## See also

Other PIPET:
[`DoPIPET()`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md)
