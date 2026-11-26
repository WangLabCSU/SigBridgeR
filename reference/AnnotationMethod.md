# AnnotationMethod: Annotation Method Class for Cell Type Identification

`AnnotationMethod` is the S7 class for cell type annotation methods. It
extends
[Method](https://wanglabcsu.github.io/sigbridger/reference/Method.md)
with Python environment support, enabling annotation methods that depend
on Python-based tools (e.g., SingleR, scArches).

## Usage

``` r
AnnotationMethod(
  method_name = character(),
  executor = function() {
     NULL
 },
  pkg_name = character(),
  method_version = NULL,
  py_env = NULL,
  sigbridger_version = get_pkg_version()
)
```

## Arguments

- method_name:

  `character`. A string identifying the method (e.g., `"Scissor"`,
  `"scPAS"`).

- executor:

  `function`. The function that performs the actual computation. This is
  the core workhorse of each method.

- pkg_name:

  The name of the R (Python) package when using this annotation method.
  When `method_version` is `NULL`, the method version is resolved from
  the R package version according to `pkg_name`.

- method_version:

  `character`. A string indicating the version of the underlying
  software package.

- py_env:

  `character` or `NULL`. A conda environment name or Python path for
  Python-based annotation tools. When `NULL` (default), the method
  version is resolved from an R package.

- sigbridger_version:

  Package Version of SigBridgeR

## Value

An `AnnotationMethod` S7 object.

## Details

This class inherits from
[Method](https://wanglabcsu.github.io/sigbridger/reference/Method.md)
and adds a `py_env` property for specifying a conda environment or
Python path. When `py_env` is `NULL`, the method version is resolved
from an R package via
[`r_pkg_version()`](https://wanglabcsu.github.io/sigbridger/reference/r_pkg_version.md).
When a Python environment is specified, the version is resolved via
[`py_pkg_version()`](https://wanglabcsu.github.io/sigbridger/reference/py_pkg_version.md).

## See also

Other S7-Classes:
[`BenchAccConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAccConfig.md),
[`BenchConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchConfig.md),
[`BenchGSEAConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEAConfig.md),
[`BenchGeneScoreConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfig.md),
[`BenchMemoryConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfig.md),
[`ScreenMethodCache()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodCache.md),
[`ScreenMethodConfig()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodConfig.md),
[`ScreenMethodResult()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodResult.md)
