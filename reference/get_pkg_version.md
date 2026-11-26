# Get SigBridgeR Package Version

Returns the installed version of the SigBridgeR package. When the
package is loaded via
[`devtools::load_all()`](https://devtools.r-lib.org/reference/load_all.html)
(i.e., in development mode), falls back to reading the `Version` field
from the `DESCRIPTION` file in the working directory.

## Usage

``` r
get_pkg_version()
```

## Value

A character string with the SigBridgeR version.

## See also

Other get-pkg-version:
[`py_pkg_version()`](https://wanglabcsu.github.io/sigbridger/reference/py_pkg_version.md),
[`r_pkg_version()`](https://wanglabcsu.github.io/sigbridger/reference/r_pkg_version.md)
