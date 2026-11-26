# Build GSEA Benchmark Controls

Collects the arguments of a GSEA benchmark into the `bench_config` list
of
[`BenchGSEAConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEAConfig.md).
Arguments matching
[`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html) formals are
handed over to it; all other named arguments are kept as-is.

## Usage

``` r
BenchGSEAConfigSet(sc_data, pathway, group = NULL, ...)
```

## Arguments

- sc_data:

  A
  [Seurat](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
  object to be screened. It defines the cells that `group` assigns
  labels to.

- pathway:

  A named `list` of gene sets, each element a character vector of gene
  symbols. Passed to the `pathways` argument of
  [`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html).

- group:

  `character` or `NULL`. Per-cell group labels, one per cell in
  `sc_data`. The vector is stored in `bench_config` and validated by
  [`BenchGSEAConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEAConfig.md).
  Default: `NULL`.

- ...:

  Additional named arguments, typically further
  [`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html) controls
  (e.g. `nperm`, `scoreType`). Case variants of the function's own
  arguments are resolved to their canonical names.

## Value

A list of parameters for
[`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html) and
[`BenchGSEAConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEAConfig.md).

## See also

Other ConfigSet:
[`BenchGeneScoreConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfigSet.md),
[`BenchMemoryConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfigSet.md),
[`MethodParamSet()`](https://wanglabcsu.github.io/sigbridger/reference/MethodParamSet.md)

## Examples

``` r
if (FALSE) { # \dontrun{
config <- BenchGSEAConfigSet(
  sc_data = sc,
  pathway = fgsea::examplePathways,
  group = sc$cell_type
)
} # }
```
