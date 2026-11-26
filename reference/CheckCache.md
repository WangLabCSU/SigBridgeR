# Check Cache Configuration Consistency

Validates whether the cached configuration matches the current
parameters. This function reads a cache metadata JSON file and compares
it with the provided screen method, phenotype class, label type, and
parameters.

## Usage

``` r
CheckCache(cache_config, ...)
```

## Arguments

- cache_config:

  A `ScreenMethodConfig` or `ScreenMethodCache` S7 object containing the
  current configuration to validate.

- ...:

  Additional arguments to be passed to the method:

  For `ScreenMethodConfig`:

  - `path`: `character`. Path to the cache directory or directly to the
    `cache_config.json` file. Required for `ScreenMethodConfig`

## Value

Invisible `TRUE` if the cache configuration is consistent with the
current parameters. Otherwise, aborts with an error message displaying
the differences.

## Details

The function compares each field in the current configuration against
the cached metadata. If any field differs (method name, method version,
phenotype class, label type, or parameters), the function aborts with a
formatted error message showing which parameter differs and both the
provided and actual values.

## Methods

`CheckCache` is an S7 generic with methods available for the following
classes:

- `ANY`

- [`SigBridgeR::ScreenMethodCache`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodCache.md)

- [`SigBridgeR::ScreenMethodConfig`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodConfig.md)

## See also

Other cache_config:
[`CacheSetHere()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSetHere.md),
[`CacheSysCall()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSysCall.md),
[`ChooseCache()`](https://wanglabcsu.github.io/sigbridger/reference/ChooseCache.md),
[`LoadCache()`](https://wanglabcsu.github.io/sigbridger/reference/LoadCache.md),
[`WriteCache()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCache.md),
[`WriteCacheMeta()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCacheMeta.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Check cache consistency for a ScreenMethodConfig
config <- ScreenMethodConfig(
  method_name = "Scissor",
  param = list(alpha = 0.05)
)
CheckCache(config, path = "cache/cache_config.json")

# Check using a ScreenMethodCache object
cache <- ScreenMethodCache(
  cache_path = "cache/",
  cache_config_path = "cache/cache_config.json",
  cache_data = list(),
  screen_method_config = config
)
CheckCache(cache)
} # }
```
