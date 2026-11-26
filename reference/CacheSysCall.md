# System-Level Cache Load and Save

High-level wrapper for loading and saving cached screening results.
Combines cache path resolution (via
[`CacheSetHere()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSetHere.md)),
configuration validation (via
[`CheckCache()`](https://wanglabcsu.github.io/sigbridger/reference/CheckCache.md)),
and data I/O (via
[`LoadCache()`](https://wanglabcsu.github.io/sigbridger/reference/LoadCache.md)
or
[`WriteCache()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCache.md))
into a single call.

## Usage

``` r
CacheSysCall(
  mode = c("load", "save"),
  path = NULL,
  cache,
  verbose = TRUE,
  timestamp = NULL,
  ...
)
```

## Arguments

- mode:

  `character`. Either `"load"` or `"save"`.

- path:

  `character` or `NULL`. Path passed to
  [`CacheSetHere()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSetHere.md).
  In save mode a `NULL` value falls back to `cache@cache_path`. In load
  mode it is mandatory: a `NULL` value aborts, no cache is looked up
  implicitly. Default: `NULL`.

- cache:

  A `ScreenMethodConfig` or `ScreenMethodCache` S7 object. For save
  mode, a `ScreenMethodCache` is required (contains the data to write).
  For load mode, either class is accepted.

- verbose:

  `logical`. Whether to print progress messages. Default: `TRUE`.

- timestamp:

  `character` or `NULL`. Optional timestamp string for the new cache
  directory in save mode. Passed to
  [`CacheSetHere()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSetHere.md).

- ...:

  Additional arguments passed to
  [`WriteCache()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCache.md)
  in save mode (e.g., `format`, `max_rows`, `max_cols`).

## Value

In `"load"` mode, a named list of loaded cache data objects. In `"save"`
mode, the cache data invisibly.

## Details

`CacheSysCall` is the primary entry point for cache operations in the
SigBridgeR pipeline. It dispatches to either `CacheSysCall.load()` or
`CacheSysCall.save()` based on the `mode` argument.

### Load mode

In `"load"` mode, the function:

1.  Resolves the cache path via
    [`CacheSetHere()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSetHere.md).

2.  Validates cache consistency via
    [`CheckCache()`](https://wanglabcsu.github.io/sigbridger/reference/CheckCache.md).

3.  Loads each cached data object via
    [`LoadCache()`](https://wanglabcsu.github.io/sigbridger/reference/LoadCache.md)
    and returns them as a named list.

### Save mode

In `"save"` mode, the function:

1.  Resolves and creates the cache directory via
    [`CacheSetHere()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSetHere.md).

2.  Writes the cache data and metadata via
    [`WriteCache()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCache.md).

3.  Returns the cache data invisibly.

## See also

Other cache_config:
[`CacheSetHere()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSetHere.md),
[`CheckCache()`](https://wanglabcsu.github.io/sigbridger/reference/CheckCache.md),
[`ChooseCache()`](https://wanglabcsu.github.io/sigbridger/reference/ChooseCache.md),
[`LoadCache()`](https://wanglabcsu.github.io/sigbridger/reference/LoadCache.md),
[`WriteCache()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCache.md),
[`WriteCacheMeta()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCacheMeta.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Load cached results
cache <- ScreenMethodCache(
  cache_path = "Scissor_res/",
  cache_config_path = "Scissor_res/cache_config.json",
  cache_data = list(),
  screen_method_config = ScreenMethodConfig(
    screen_method = "Scissor",
    param = list(alpha = 0.05)
  )
)
result <- CacheSysCall(
  "load",
  path = "Scissor_res/survival_2025_01_01",
  cache = cache
)

# Save results to cache
cache <- ScreenMethodCache(
  cache_path = "Scissor_res/",
  cache_config_path = "Scissor_res/cache_config.json",
  cache_data = list(result = my_matrix),
  screen_method_config = ScreenMethodConfig(
    screen_method = "Scissor",
    param = list(alpha = 0.05)
  )
)
CacheSysCall("save", cache = cache)
} # }
```
