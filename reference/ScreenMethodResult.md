# ScreenMethodResult: Screening Method Result Container

`ScreenMethodResult` stores the results of a phenotype-associated cell
screening analysis. It wraps a Seurat object containing the screened
single-cell data and supports attaching additional dynamic properties
(e.g., screening statistics, metadata) via `...` in the constructor.

## Usage

``` r
ScreenMethodResult(scRNA_data, ..., sigbridger_version = get_pkg_version())
```

## Arguments

- scRNA_data:

  A
  [Seurat](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
  object containing the screening results. This is the core required
  property of the class.

- ...:

  Named arguments to set as additional S7 properties on the result
  object (e.g., screening statistics, method-specific outputs). Stored
  as a list in the `extra` property.

- sigbridger_version:

  Package Version of SigBridgeR

## Value

A `ScreenMethodResult` S7 object.

## Details

This class is defined using the S7 object-oriented system and inherits
from
[SigBridgeRBase](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeRBase.md).
The core required property is `scRNA_data`, which must be a
[Seurat](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
object containing the screening results.

Additional properties can be attached dynamically via the constructor's
`...` argument. Named arguments are set as S7 properties on the returned
object, allowing flexible result storage without subclassing.

### Validator

The class validator enforces:

- `scRNA_data` must be present on the object.

- `scRNA_data` must be a
  [Seurat](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
  object.

## See also

Other S7-Classes:
[`AnnotationMethod()`](https://wanglabcsu.github.io/sigbridger/reference/AnnotationMethod.md),
[`BenchAccConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAccConfig.md),
[`BenchConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchConfig.md),
[`BenchGSEAConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEAConfig.md),
[`BenchGeneScoreConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfig.md),
[`BenchMemoryConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfig.md),
[`ScreenMethodCache()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodCache.md),
[`ScreenMethodConfig()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodConfig.md)
