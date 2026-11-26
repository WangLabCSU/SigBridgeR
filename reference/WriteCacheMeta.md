# Write Cache Metadata

Writes metadata for SigBridgeR cache identification to a JSON file. The
metadata includes general system information, cache configuration, and
an optional user-provided description.

## Usage

``` r
WriteCacheMeta(cache_config, ...)
```

## Arguments

- cache_config:

  A `ScreenMethodConfig` or `ScreenMethodCache` S7 object containing the
  cache configuration to write.

- ...:

  Additional arguments passed to methods.

  For `ScreenMethodConfig`:

  - `file`: `character(1)`. Path to the JSON metadata file to be
    written.

  - `additional_description`: `character` or `NULL`. An optional string
    describing the cache entry. Default: `NULL`.

  - `verbose`: `logical`. Whether to print progress messages. Default:
    `TRUE`.

  For `ScreenMethodCache`:

  - `additional_description`: `character` or `NULL`. An optional string
    describing the cache entry. Default: `NULL`.

  - `verbose`: `logical`. Whether to print progress messages. Default:
    `TRUE`.

## Value

Invisible. Returns the metadata list that was written to the file.

## Details

This function creates a JSON metadata file containing:

- **`general`**: File path, OS information, timestamp, R version, and
  SigBridgeR version

- **`config`**: Screen method, phenotype class, label type, and
  parameters

- **`description`**: Additional user-provided description

The JSON file includes a header comment indicating it is required for
cache identification and should not be modified.

Only for SigBridgeR version \>= 3.7.0.

## Methods

`WriteCacheMeta` is an S7 generic with methods available for the
following classes:

- `ANY`

- [`SigBridgeR::ScreenMethodCache`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodCache.md)

- [`SigBridgeR::ScreenMethodConfig`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodConfig.md)

## See also

Other cache_config:
[`CacheSetHere()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSetHere.md),
[`CacheSysCall()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSysCall.md),
[`CheckCache()`](https://wanglabcsu.github.io/sigbridger/reference/CheckCache.md),
[`ChooseCache()`](https://wanglabcsu.github.io/sigbridger/reference/ChooseCache.md),
[`LoadCache()`](https://wanglabcsu.github.io/sigbridger/reference/LoadCache.md),
[`WriteCache()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCache.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Write metadata for a ScreenMethodConfig
config <- ScreenMethodConfig(
  method_name = "Scissor",
  param = list(alpha = 0.05)
)
WriteCacheMeta(config, file = "cache/cache_config.json")

# With an additional description
WriteCacheMeta(config,
  file = "cache/cache_config.json",
  additional_description = "Scissor screening for survival phenotype"
)
} # }
```
