# Get Installed Version of a Python Package in a Conda Environment

Queries the installed version of a Python package in a specified conda
environment by running `conda list -n <env> <pkg>`. Returns `NULL` if
the package is not found or conda is unavailable.

## Usage

``` r
py_pkg_version(pkg_name, conda_env)
```

## Arguments

- pkg_name:

  A character string. The name of the Python package.

- conda_env:

  A character string. The name of the conda environment.

## Value

A character string with the version, or `NULL` if unavailable.

## See also

Other get-pkg-version:
[`get_pkg_version()`](https://wanglabcsu.github.io/sigbridger/reference/get_pkg_version.md),
[`r_pkg_version()`](https://wanglabcsu.github.io/sigbridger/reference/r_pkg_version.md)
