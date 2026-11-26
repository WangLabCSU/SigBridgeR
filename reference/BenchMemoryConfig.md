# BenchMemoryConfig: Time and Memory Benchmark Configuration

`BenchMemoryConfig` stores the configuration of a time and memory
benchmark run performed by
[`BenchMemory()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemory.md).
It extends
[BenchConfig](https://wanglabcsu.github.io/sigbridger/reference/BenchConfig.md)
with the benchmarking controls collected in `bench_config`; the inputs
shared by every screening method (`matched_bulk`, `sc_data`,
`phenotype`, `label_type`, `phenotype_class`, `screen_method`) and the
method specific arguments collected in `method_param` are inherited from
[BenchConfig](https://wanglabcsu.github.io/sigbridger/reference/BenchConfig.md).

## Usage

``` r
BenchMemoryConfig(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = "{method_name}",
  phenotype_class = c("binary", "survival", "continuous"),
  screen_method = names(ScreenStrategy),
  method_param = MethodParamSet(),
  bench_config = BenchMemoryConfigSet(),
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

  Named `list` of benchmarking controls (rounds, iterations, order,
  budgets, ...) consumed by
  [`BenchMemory()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemory.md).
  Default: [`list()`](https://rdrr.io/r/base/list.html), meaning all
  defaults. Unknown fields abort.

- ...:

  Must be empty. Method specific arguments belong in `method_param`;
  unlike the other benchmark configurations nothing is stored in
  `param`.

## Value

A `BenchMemoryConfig` S7 object.

## See also

Other S7-Classes:
[`AnnotationMethod()`](https://wanglabcsu.github.io/sigbridger/reference/AnnotationMethod.md),
[`BenchAccConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAccConfig.md),
[`BenchConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchConfig.md),
[`BenchGSEAConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEAConfig.md),
[`BenchGeneScoreConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfig.md),
[`ScreenMethodCache()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodCache.md),
[`ScreenMethodConfig()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodConfig.md),
[`ScreenMethodResult()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodResult.md)
