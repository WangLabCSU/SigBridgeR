# Resolve and Set Up Cache Directory

Resolves a user-specified path to the cache layer directory. In `"save"`
mode, the necessary directory structure is created automatically. In
`"load"` mode, an existing cache directory is identified.

## Usage

``` r
CacheSetHere(
  path,
  cache_config,
  mode = c("load", "save"),
  timestamp = NULL,
  ...
)
```

## Arguments

- path:

  `character`. User-specified path to a cache, root, or parent
  directory.

- cache_config:

  A `ScreenMethodConfig` or `ScreenMethodCache` S7 object. `method_name`
  and `phenotype_class` are extracted from this object to construct the
  cache directory path.

- mode:

  `character`. Either `"load"` or `"save"`. Default: `"load"`.

- timestamp:

  `character` or `NULL`. Optional timestamp string for the new cache
  directory in save mode. When `NULL`, defaults to
  `format(Sys.time(), "%Y%m%d%H%M%S")`.

- ...:

  Additional arguments (must be empty, checked by
  [`rlang::check_dots_empty0()`](https://rlang.r-lib.org/reference/check_dots_empty0.html)).

## Value

A `character` string: the absolute path to the cache directory.

## Details

The function detects which cache layer the given path belongs to:

- **Cache layer**: The path directly contains a `cache_config.json`
  file. In save mode this triggers an error (recursive caching is not
  supported); in load mode the path is returned as-is.

- **Root layer**: The path's basename matches `{method_name}_res` (e.g.,
  `"Scissor_res"`). In save mode a new cache subdirectory is created; in
  load mode
  [`ChooseCache()`](https://wanglabcsu.github.io/sigbridger/reference/ChooseCache.md)
  is called to select an existing cache.

- **Parent layer**: The path is a parent directory. A
  `{method_name}_res` subdirectory is created (save) or located (load)
  underneath it.

In save mode, the cache directory is named
`{phenotype_class}_{timestamp}` (e.g., `"survival_20250730120000"`).

## See also

Other cache_config:
[`CacheSysCall()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSysCall.md),
[`CheckCache()`](https://wanglabcsu.github.io/sigbridger/reference/CheckCache.md),
[`ChooseCache()`](https://wanglabcsu.github.io/sigbridger/reference/ChooseCache.md),
[`LoadCache()`](https://wanglabcsu.github.io/sigbridger/reference/LoadCache.md),
[`WriteCache()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCache.md),
[`WriteCacheMeta()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCacheMeta.md)

## Examples

``` r
if (FALSE) { # \dontrun{
config <- ScreenMethodConfig(
  method_name = "Scissor",
  param = list(alpha = 0.05)
)

# Save: create a new cache under the root
CacheSetHere("Scissor_res", config, mode = "save")

# Load: select an existing cache from the root
CacheSetHere("Scissor_res", config, mode = "load")

# Load: point directly to a specific cache
CacheSetHere("Scissor_res/survival_207701011212", config, mode = "load")
} # }
```
