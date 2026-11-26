# Choose a Cache Directory

Selects a cache directory from the specified path. If only one cache
subdirectory exists, it is returned automatically; if multiple are
found, an interactive menu is presented for the user to choose.

## Usage

``` r
ChooseCache(directory)
```

## Arguments

- directory:

  `character`. Path to the parent directory containing cache
  subdirectories. The directory name must match `{method_name}_res` for
  a registered screening method.

## Value

A `character` string: the path to the selected cache directory.

## Details

The function expects the parent directory name to follow the convention
`{method_name}_res` (e.g., `"Scissor_res"`, `"scPAS_res"`). It validates
the directory name against the registered screening methods in
`ScreenStrategy`.

If no cache subdirectories are found, the function aborts with an error.
If exactly one subdirectory exists, it is returned without prompting. If
multiple subdirectories exist,
[`utils::menu()`](https://rdrr.io/r/utils/menu.html) is used to present
an interactive selection menu.

## See also

Other cache_config:
[`CacheSetHere()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSetHere.md),
[`CacheSysCall()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSysCall.md),
[`CheckCache()`](https://wanglabcsu.github.io/sigbridger/reference/CheckCache.md),
[`LoadCache()`](https://wanglabcsu.github.io/sigbridger/reference/LoadCache.md),
[`WriteCache()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCache.md),
[`WriteCacheMeta()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCacheMeta.md)

## Examples

``` r
if (FALSE) { # \dontrun{
cache_path <- ChooseCache("./Scissor_res")
} # }
```
