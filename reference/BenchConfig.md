# BenchConfig: Base Configuration for Benchmarks

`BenchConfig` is the common parent of the benchmark configurations in
SigBridgeR. It stores the inputs shared by every screening method
(`matched_bulk`, `sc_data`, `phenotype`, `label_type`,
`phenotype_class`, `screen_method`), the method specific arguments
collected in `method_param`, and an optional `param` list used by
subclasses that forward extra arguments.

## Usage

``` r
BenchConfig(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = "{method_name}",
  phenotype_class = c("binary", "survival", "continuous"),
  screen_method = names(ScreenStrategy),
  method_param = MethodParamSet(),
  sigbridger_version = get_pkg_version(),
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

- sigbridger_version:

  Package Version of SigBridgeR

- ...:

  unused

## Value

A `BenchConfig` S7 object.

## Details

A single configuration may describe several screening runs: the
benchmark screens the same data with each method listed in
`screen_method`, one method at a time.

- `screen_method` may only contain methods registered in
  [ScreenStrategy](https://wanglabcsu.github.io/sigbridger/reference/ScreenStrategy.md)
  (e.g. `"Scissor"`, `"scPAS"`, `"LP_SGL"`); the constructor of the
  concrete configuration aborts on any unknown name. It also aborts when
  the selected methods have no phenotype class in common.

- `method_param` is a named list with one entry per method, named
  `"<method>_param"` (e.g. `list(Scissor_param = list(alpha = 0.05))`).
  Each entry is spliced into the arguments of that method's executor,
  after its
  [ScreenMethod](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethod.md)
  mapper has run. Use
  [`MethodParamSet()`](https://wanglabcsu.github.io/sigbridger/reference/MethodParamSet.md)
  to build the full set with empty defaults, so unspecified methods
  simply run with their own defaults.

- `label_type` is a length-one template. The default `"{method_name}"`
  is resolved by the benchmark with the name of the method currently
  being benchmarked, so a single configuration can serve methods whose
  executors expect differently named phenotype labels.

## See also

Other S7-Classes:
[`AnnotationMethod()`](https://wanglabcsu.github.io/sigbridger/reference/AnnotationMethod.md),
[`BenchAccConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAccConfig.md),
[`BenchGSEAConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEAConfig.md),
[`BenchGeneScoreConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfig.md),
[`BenchMemoryConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfig.md),
[`ScreenMethodCache()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodCache.md),
[`ScreenMethodConfig()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodConfig.md),
[`ScreenMethodResult()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodResult.md)
