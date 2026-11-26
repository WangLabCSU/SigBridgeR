# Python Module Registry

An internal environment holding the Python modules bundled under the
package's `inst/python/` directory. Modules are loaded lazily at package
startup via
[`reticulate::import_from_path()`](https://rstudio.github.io/reticulate/reference/import.html)
and can be accessed as `PyModule$<module>`.

## Usage

``` r
PyModule
```

## Format

An environment.

## Details

The following Python modules are currently bundled with the package:

- `celltypist_annotate`

These modules provide the Python backend used by R functions such as
[`CellTypistAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/CellTypistAnnotate.md).
