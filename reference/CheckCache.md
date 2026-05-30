# Check Cache Configuration Consistency

Validates whether the cached configuration matches the current
parameters. This function reads a cache metadata JSON file and compares
it with the provided screen method, phenotype class, label type, and
parameters.

## Usage

``` r
CheckCache(
  path,
  screen_method,
  phenotype_class = c("binary", "survival", "continuous"),
  label_type = screen_method,
  params,
  ...
)
```

## Arguments

- path:

  Character string specifying the path to the cache directory or
  directly to the `cache_config.json` file.

- screen_method:

  Character string indicating the screening method used.

- phenotype_class:

  Character vector specifying the phenotype class type. Must be one of
  `"binary"`, `"survival"`, or `"continuous"`.

- label_type:

  Character string specifying the label type. Defaults to the value of
  `screen_method`.

- params:

  List of parameters used for the screening method.

- ...:

  Additional arguments (must be empty; raises error if provided).

## Value

Returns `invisible(TRUE)` if the cache configuration is consistent with
the current parameters. Otherwise, aborts with an error message
displaying the differences.

## See also

Other cache_config:
[`CacheSetHere()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSetHere.md),
[`ChooseCache()`](https://wanglabcsu.github.io/sigbridger/reference/ChooseCache.md),
[`LoadCache()`](https://wanglabcsu.github.io/sigbridger/reference/LoadCache.md),
[`WriteCache()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCache.md),
[`WriteCacheMeta()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCacheMeta.md)
