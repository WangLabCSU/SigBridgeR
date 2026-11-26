# Write an R Object to Cache

Writes an R object to the cache directory. Uses `qs2` for fast
serialization when available, falling back to `rdata` for general
objects. Small data frames and data.tables are saved as CSV via
[`data.table::fwrite()`](https://rdrr.io/pkg/data.table/man/fwrite.html)
for easy inspection.

## Usage

``` r
WriteCache(cache, ...)
```

## Arguments

- cache:

  A `ScreenMethodConfig` or `ScreenMethodCache` S7 object containing the
  data and configuration to write.

- ...:

  Additional arguments passed to methods.

  For `ScreenMethodConfig`:

  - `file`: `character(1)`. Path to the output `cache_config.json` file.

  - `additional_description`: `character` or `NULL`. An optional string
    describing the cache entry. Default: `NULL`.

  - `verbose`: `logical`. Whether to print progress messages. Default:
    `TRUE`.

  For `ScreenMethodCache`:

  - `format`: `character`. Cache format. One of `"auto"`, `"qs2"`,
    `"qdata"`, `"rdata"`, or `"csv"`. When `"auto"` default, the format
    is detected per object.

  - `max_rows`: `integer`. Maximum number of rows to consider a
    data.frame small for CSV output. Default: `1000L`.

  - `max_cols`: `integer`. Maximum number of columns to consider a
    data.frame small for CSV output. Default: `20L`.

  - `additional_description`: `character` or `NULL`. An optional string
    describing the cache entry. Default: `NULL`.

  - `verbose`: `logical`. Whether to print progress messages. Default:
    `TRUE`.

## Value

Invisible. Returns the absolute path(s) to the written cache file(s).

## Details

The `ScreenMethodCache` method writes both the cache metadata (via
[`WriteCacheMeta()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCacheMeta.md))
and the cached data objects to disk. The `ScreenMethodConfig` method
only writes the metadata.

### Format detection (ScreenMethodCache)

When `format = "auto"` (default), the format for each data object is
detected automatically:

- If the object is a 2D structure (data.frame/matrix) with at most
  `max_rows` rows and `max_cols` columns, CSV format is used.

- If `qs2` is installed and the object is a list, vector, data.frame, or
  matrix, `qdata` format is used.

- If `qs2` is installed but the object is not qdata-compatible, `qs2`
  format is used.

- Otherwise, `rdata` format is used as a fallback.

## Methods

`WriteCache` is an S7 generic with methods available for the following
classes:

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
[`WriteCacheMeta()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCacheMeta.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Write using a ScreenMethodCache
cache <- ScreenMethodCache(
  cache_path = "Scissor_res/survival_2025_01_01/",
  cache_config_path = "Scissor_res/survival_2025_01_01/cache_config.json",
  cache_data = list(result = my_matrix),
  screen_method_config = config
)
WriteCache(cache)

# Force CSV format for small data frames
WriteCache(cache, format = "csv")
} # }
```
