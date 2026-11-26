# ScreenMethodConfig: Screening Method Configuration

`ScreenMethodConfig` stores the configuration for a phenotype-associated
cell screening method. It is used both to define screening parameters
and as a cache result descriptor.

## Usage

``` r
ScreenMethodConfig(
  method_name,
  param,
  method_version = NULL,
  phenotype_class = c("binary", "continuous", "survival"),
  label_type = character(0L),
  sigbridger_version = get_pkg_version()
)
```

## Arguments

- method_name:

  `character`. The name of the screening method (e.g., `"Scissor"`,
  `"scPAS"`). Must match a key in
  [ScreenStrategy](https://wanglabcsu.github.io/sigbridger/reference/ScreenStrategy.md).

- param:

  `list`. A named list of parameters passed to the screening function.
  Must be non-empty.

- method_version:

  `character` or `NULL`. The version of the underlying software package.
  When `NULL`, automatically resolved via `r_pkg_version(method_name)`.

- phenotype_class:

  `character`. The phenotype classification type(s). One or more of
  `"binary"`, `"continuous"`, `"survival"`. Default: all three.

- label_type:

  `character`. The label type for the screening method. Default:
  `character(0L)`.

- sigbridger_version:

  Package Version of SigBridgeR

## Value

A `ScreenMethodConfig` S7 object.

## Details

This class is defined using the S7 object-oriented system and inherits
from
[SigBridgeRBase](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeRBase.md).
It validates that `method_name` is a registered screening method in
[ScreenStrategy](https://wanglabcsu.github.io/sigbridger/reference/ScreenStrategy.md)
and that `phenotype_class` is one of the supported types.

### Validator

The class validator enforces:

- `param` must be a non-empty named list.

- `phenotype_class` must be one or more of `"binary"`, `"continuous"`,
  `"survival"`.

- `method_name` must match a key in
  [ScreenStrategy](https://wanglabcsu.github.io/sigbridger/reference/ScreenStrategy.md).

## See also

Other S7-Classes:
[`AnnotationMethod()`](https://wanglabcsu.github.io/sigbridger/reference/AnnotationMethod.md),
[`BenchAccConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAccConfig.md),
[`BenchConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchConfig.md),
[`BenchGSEAConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEAConfig.md),
[`BenchGeneScoreConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfig.md),
[`BenchMemoryConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfig.md),
[`ScreenMethodCache()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodCache.md),
[`ScreenMethodResult()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodResult.md)
