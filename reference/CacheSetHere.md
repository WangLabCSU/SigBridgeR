# Resolve and Set Up Cache Directory

Resolves a user-specified path to the cache layer directory. In `"save"`
mode, the necessary directory structure is created automatically. In
`"load"` mode, an existing cache directory is identified.

## Usage

``` r
CacheSetHere(
  path,
  screen_method,
  phenotype_class = c("binary", "survival", "continuous"),
  mode = c("load", "save"),
  timestamp = NULL,
  ...
)
```

## Arguments

- path:

  Character string. User-specified path.

- screen_method:

  Character string. Screening method name. Must match a key in
  [`ScreenStrategy`](https://wanglabcsu.github.io/sigbridger/reference/ScreenStrategy.md).

- phenotype_class:

  Character string. Phenotype class. One of `"binary"`, `"survival"`, or
  `"continuous"`.

- mode:

  Character string. Either `"load"` (default) or `"save"`.

- timestamp:

  Character string. Optional timestamp string for the new cache
  directory in save mode. Defaults to
  `format(Sys.time(), "%Y_%m_%d_%H%M")`.

- ...:

  Additional arguments (must be empty).

## Value

Character string. The absolute path to the cache layer directory.

## See also

Other cache_config:
[`CheckCache()`](https://wanglabcsu.github.io/sigbridger/reference/CheckCache.md),
[`ChooseCache()`](https://wanglabcsu.github.io/sigbridger/reference/ChooseCache.md),
[`LoadCache()`](https://wanglabcsu.github.io/sigbridger/reference/LoadCache.md),
[`WriteCache()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCache.md),
[`WriteCacheMeta()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCacheMeta.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Save: create a new cache under the root
CacheSetHere("Scissor_res", "Scissor", "survival", mode = "save")

# Load: select an existing cache from the root
CacheSetHere("Scissor_res", "Scissor", "survival", mode = "load")

# Load: point directly to a specific cache
CacheSetHere("Scissor_res/survival_202512011212", "Scissor", "survival", mode = "load")
} # }
```
