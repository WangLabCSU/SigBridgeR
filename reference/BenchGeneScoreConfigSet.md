# Build Gene Set Scoring Benchmark Controls

Collects the arguments of a gene set scoring benchmark into the
`bench_config` list of
[`BenchGeneScoreConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfig.md).
`method` selects the scoring backend, and additional named arguments are
forwarded to that backend when
[`BenchGeneScore()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScore.md)
computes the scores.

## Usage

``` r
BenchGeneScoreConfigSet(
  method = c("ssGSEA", "AUCell", "UCell", "AddModuleScore"),
  ...
)
```

## Arguments

- method:

  Character. Scoring backend to benchmark, one of `"ssGSEA"`,
  `"AUCell"`, `"UCell"`, `"AddModuleScore"`; The package backing the
  selected backend is checked with
  [`rlang::check_installed()`](https://rlang.r-lib.org/reference/is_installed.html).
  Default: `"ssGSEA"`.

- ...:

  Additional named arguments forwarded to the selected scoring backend.
  They are matched against the backend's formal arguments and, when a
  method needs two calls, split over its stages (e.g.
  [`GSVA::ssgseaParam()`](https://rdrr.io/pkg/GSVA/man/ssgseaParam-class.html)
  and [`GSVA::gsva()`](https://rdrr.io/pkg/GSVA/man/gsva.html) for
  `"ssGSEA"`). Case variants of the function's own arguments are
  resolved to their canonical names.

## Value

A named `list` holding `method` and the extra backend arguments, to be
used as the `bench_config` of
[`BenchGeneScoreConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfig.md).

## See also

[`BenchGeneScoreConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfig.md),
[`BenchGeneScore()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScore.md)

Other ConfigSet:
[`BenchGSEAConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEAConfigSet.md),
[`BenchMemoryConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfigSet.md),
[`MethodParamSet()`](https://wanglabcsu.github.io/sigbridger/reference/MethodParamSet.md)

## Examples

``` r
BenchGeneScoreConfigSet(method = "ssGSEA")
#> $method
#> [1] "ssGSEA"
#> 
```
