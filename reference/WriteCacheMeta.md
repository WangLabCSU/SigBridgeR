# Write Cache Metadata

Writes metadata for SigBridgeR cache identification to a JSON file.

## Usage

``` r
WriteCacheMeta(
  file = "cache_config.json",
  screen_method,
  phenotype_class = c("binary", "survival", "continuous"),
  label_type = screen_method,
  params,
  additional_description = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- file:

  Character string. Path to the JSON metadata file to be written.

- screen_method:

  Character string. Screening method name. Must match one of the keys in
  `ScreenStrategy`.

- phenotype_class:

  Character string. Type of phenotype data. One of `"binary"`,
  `"survival"`, or `"continuous"`.

- label_type:

  Character string. Type of label for the screening. Defaults to
  `screen_method`.

- params:

  List. Parameters used for the screening configuration.

- additional_description:

  Character string. Optional additional description text to include in
  the metadata. Default is `NULL`.

- verbose:

  Logical. Whether to print metadata preview for debugging. Default is
  `FALSE`.

- ...:

  Additional arguments (must be empty, checked by
  [`rlang::check_dots_empty0()`](https://rlang.r-lib.org/reference/check_dots_empty0.html)).

## Value

Invisible. Returns the metadata list that was written to the file.

## Details

This function creates a JSON metadata file containing:

- general:

  File path, OS information, timestamp, R version, and SigBridgeR
  version

- config:

  Screen method, phenotype class, label type, and parameters

- description:

  Additional user-provided description

Only for SigBridgeR version \>= 3.7.0

The JSON file includes a header comment indicating it is required for
cache identification and should not be modified.

## See also

Other cache_config:
[`CacheSetHere()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSetHere.md),
[`CheckCache()`](https://wanglabcsu.github.io/sigbridger/reference/CheckCache.md),
[`ChooseCache()`](https://wanglabcsu.github.io/sigbridger/reference/ChooseCache.md),
[`LoadCache()`](https://wanglabcsu.github.io/sigbridger/reference/LoadCache.md),
[`WriteCache()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCache.md)

## Examples

``` r
if (FALSE) { # \dontrun{
WriteCacheMeta(
  file = "cache/meta.json",
  screen_method = "DEGAS",
  phenotype_class = "binary",
  params = list(n_top = 100),
  additional_description = "Test run"
)
} # }
```
