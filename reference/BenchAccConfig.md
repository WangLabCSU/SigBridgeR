# BenchAccConfig: Accuracy Benchmark Configuration

`BenchAccConfig` stores the configuration of an accuracy benchmark run
performed by
[`BenchAcc()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAcc.md).
It extends
[BenchConfig](https://wanglabcsu.github.io/sigbridger/reference/BenchConfig.md)
with the ground-truth cell labels collected in `bench_config`; the
inputs shared by every screening method (`matched_bulk`, `sc_data`,
`phenotype`, `label_type`, `phenotype_class`, `screen_method`) and the
method specific arguments collected in `method_param` are inherited from
[BenchConfig](https://wanglabcsu.github.io/sigbridger/reference/BenchConfig.md).

## Usage

``` r
BenchAccConfig(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = "{method_name}",
  phenotype_class = c("binary", "survival", "continuous"),
  screen_method = names(ScreenStrategy),
  method_param = MethodParamSet(),
  bench_config = character(),
  ...
)
```

## Arguments

- matched_bulk:

  `matrix` or `data.frame`. Bulk expression data (genes x samples)
  matching `phenotype`. Shared by all methods.

- sc_data:

  A
  [Seurat](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
  object to be screened.

- phenotype:

  Phenotype data, either a named vector or a patient survival data
  frame. See
  [`Screen()`](https://wanglabcsu.github.io/sigbridger/reference/Screen.md).

- label_type:

  `character` or `NULL`. Length-one glue-style template for the
  phenotype label type, with `{method_name}` available. When `NULL` or
  not of length one, it falls back to `"{method_name}"` with an
  informational message.

- phenotype_class:

  `character`. One or more of `"binary"`, `"continuous"`, `"survival"`.
  Default: all three.

- screen_method:

  `character`. One or more method names registered in
  [ScreenStrategy](https://wanglabcsu.github.io/sigbridger/reference/ScreenStrategy.md),
  giving the methods to benchmark and their order. Default: `"Scissor"`.

- method_param:

  Named `list` of method-specific arguments, each passed only to the
  executor of the method it is named after (e.g.
  `list(Scissor_param = list(alpha = 0.05))`). Default:
  [`MethodParamSet()`](https://wanglabcsu.github.io/sigbridger/reference/MethodParamSet.md),
  i.e. one empty entry per registered method.

- bench_config:

  `character`. Ground-truth labels used to score the screening results:
  either a length-one column name in `sc_data@meta.data`, or a per-cell
  label vector of length `ncol(sc_data)`. Must contain a `"Positive"`
  label. Default:
  [`character()`](https://rdrr.io/r/base/character.html), which is
  rejected by the validator.

- ...:

  Must be empty. Method specific arguments belong in `method_param`.

## Value

A `BenchAccConfig` S7 object.

## Details

`bench_config` supplies the reference ("ground-truth") labels that the
screening results are scored against. It is a `character` vector
interpreted in one of two ways:

- Length one: the name of a column in `sc_data@meta.data` holding the
  per-cell labels. The column must exist and contain a `"Positive"`
  level.

- Length greater than one: a per-cell label vector, one label per cell.
  It must contain a `"Positive"` level and have the same length as the
  number of cells in `sc_data`.

In both cases a `"Positive"` label is required, since
[`BenchAcc()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAcc.md)
computes binary classification metrics (TPR, FPR, TNR, FNR, F1,
Accuracy) treating `"Positive"` as the positive class. An empty
`bench_config` is rejected by the validator.

## See also

Other S7-Classes:
[`AnnotationMethod()`](https://wanglabcsu.github.io/sigbridger/reference/AnnotationMethod.md),
[`BenchConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchConfig.md),
[`BenchGSEAConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEAConfig.md),
[`BenchGeneScoreConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfig.md),
[`BenchMemoryConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfig.md),
[`ScreenMethodCache()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodCache.md),
[`ScreenMethodConfig()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodConfig.md),
[`ScreenMethodResult()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodResult.md)
