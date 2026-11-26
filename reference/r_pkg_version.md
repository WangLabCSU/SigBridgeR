# Get Installed Version of an R Package

Queries the installed version of an R package via
[`utils::packageVersion()`](https://rdrr.io/r/utils/packageDescription.html).
Returns `NULL` on error (e.g., package not installed).

## Usage

``` r
r_pkg_version(pkg_name)
```

## Arguments

- pkg_name:

  A character string. The name of the R package.

## Value

A character string with the version, or `NULL` if unavailable.

## See also

Other get-pkg-version:
[`get_pkg_version()`](https://wanglabcsu.github.io/sigbridger/reference/get_pkg_version.md),
[`py_pkg_version()`](https://wanglabcsu.github.io/sigbridger/reference/py_pkg_version.md)
