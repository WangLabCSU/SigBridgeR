# Load a Cached R Object

Reads a cached R object from disk. The file format is determined by the
file extension:

- `.qs2`: loaded via
  [`qs2::qs_read()`](https://rdrr.io/pkg/qs2/man/qs_read.html)

- `.RData` / `.rds`: loaded via
  [`readRDS()`](https://rdrr.io/r/base/readRDS.html)

- `.csv`: loaded via
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html)

If the file was originally saved in `.qs2` format but `qs2` is not
installed on the current machine, a clear error message is shown with
installation instructions.

## Usage

``` r
LoadCache(file, ...)
```

## Arguments

- file:

  Character string. Path to the cache file to load.

- ...:

  Additional arguments (must be empty, checked by
  [`rlang::check_dots_empty0()`](https://rlang.r-lib.org/reference/check_dots_empty0.html)).

## Value

The cached R object.

## See also

Other cache_config:
[`CacheSetHere()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSetHere.md),
[`CheckCache()`](https://wanglabcsu.github.io/sigbridger/reference/CheckCache.md),
[`ChooseCache()`](https://wanglabcsu.github.io/sigbridger/reference/ChooseCache.md),
[`WriteCache()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCache.md),
[`WriteCacheMeta()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCacheMeta.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- LoadCache("Scissor_res/survival_2025_01_01/result.qs2")
feature_df <- LoadCache("Scissor_res/survival_2025_01_01/feature_table.csv")
} # }
```
