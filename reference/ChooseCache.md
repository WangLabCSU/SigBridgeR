# Choose a Cache Directory

Selects a cache directory from the specified path. If only one cache
subdirectory exists, it is returned automatically; if multiple are
found, an interactive menu is presented for the user to choose.

## Usage

``` r
ChooseCache(directory, ...)
```

## Arguments

- directory:

  Path to the parent directory containing cache subdirectories.

- ...:

  Unused, reserved for future extensions.

## Value

The path to the selected cache directory (character string).

## See also

Other cache_config:
[`CheckCache()`](https://wanglabcsu.github.io/sigbridger/reference/CheckCache.md),
[`WriteCacheMeta()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCacheMeta.md)

## Examples

``` r
if (FALSE) { # interactive()
if (FALSE) { # \dontrun{
cache_path <- ChooseCache("./Scissor_res")
} # }
}
```
