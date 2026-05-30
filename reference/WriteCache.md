# Write an R Object to Cache

Writes an R object to the cache directory. Uses `qs2` for fast
serialization when available, falling back to `RDS` for general objects.
Small data frames and data.tables are saved as CSV via
[`data.table::fwrite()`](https://rdrr.io/pkg/data.table/man/fwrite.html)
for easy inspection.

## Usage

``` r
WriteCache(
  x,
  file,
  format = c("auto", "qs2", "rds", "csv"),
  max_rows = 1000L,
  max_cols = 20L,
  ...
)
```

## Arguments

- x:

  An R object to be cached.

- file:

  Character string. Path to the cache file (with or without extension).
  If no extension is given, it is appended automatically.

- format:

  Character string. Cache format. One of `"auto"`, `"qs2"`, `"rds"`, or
  `"csv"`. When `"auto"` (default):

  - If `x` is a data.frame with at most `max_rows` rows and `max_cols`
    columns, CSV format is used.

  - Otherwise, `qs2` is preferred when the package is installed, falling
    back to `RDS`.

- max_rows:

  Integer. Maximum number of rows to consider a data.frame "small".
  Default is `1000L`.

- max_cols:

  Integer. Maximum number of columns to consider a data.frame "small".
  Default is `20L`.

- ...:

  Additional arguments (must be empty, checked by
  [`rlang::check_dots_empty0()`](https://rlang.r-lib.org/reference/check_dots_empty0.html)).

## Value

Invisible. Returns the absolute path to the written cache file.

## See also

Other cache_config:
[`CacheSetHere()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSetHere.md),
[`CheckCache()`](https://wanglabcsu.github.io/sigbridger/reference/CheckCache.md),
[`ChooseCache()`](https://wanglabcsu.github.io/sigbridger/reference/ChooseCache.md),
[`LoadCache()`](https://wanglabcsu.github.io/sigbridger/reference/LoadCache.md),
[`WriteCacheMeta()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCacheMeta.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Write a matrix as qs2
WriteCache(cache_res, "Scissor_res/survival_2025_01_01/result.qs2")

# Write a small data.frame as CSV
WriteCache(small_df, "Scissor_res/survival_2025_01_01/feature_table.csv")
} # }
```
