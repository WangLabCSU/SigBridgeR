# A List of Method Parameters

Collects the method-specific arguments of a benchmark configuration.
Each method registered in
[ScreenStrategy](https://wanglabcsu.github.io/sigbridger/reference/ScreenStrategy.md)
has a corresponding argument named `"<method>_param"`, which defaults to
[`list()`](https://rdrr.io/r/base/list.html).

## Usage

``` r
MethodParamSet(
  Scissor_param = list(),
  scAB_param = list(),
  scPAS_param = list(),
  scPP_param = list(),
  DEGAS_param = list(),
  DEGASv2_param = list(),
  LP_SGL_param = list(),
  PIPET_param = list(),
  SCIPAC_param = list(),
  SIDISH_param = list(),
  TiRank_param = list(),
  PhenoMapR_param = list(),
  ...
)
```

## Arguments

- Scissor_param:

  Extra arguments for method Scissor

- scAB_param:

  Extra arguments for method scAB

- scPAS_param:

  Extra arguments for method scPAS

- scPP_param:

  Extra arguments for method scPP

- DEGAS_param:

  Extra arguments for method DEGAS

- DEGASv2_param:

  Extra arguments for method DEGASv2

- LP_SGL_param:

  Extra arguments for method LP_SGL

- PIPET_param:

  Extra arguments for method PIPET

- SCIPAC_param:

  Extra arguments for method SCIPAC

- SIDISH_param:

  Extra arguments for method SIDISH

- TiRank_param:

  Extra arguments for method TiRank

- PhenoMapR_param:

  Extra arguments for method PhenoMapR

- ...:

  Currently unused

## Value

A modified list of method-specific parameters.

## See also

Other ConfigSet:
[`BenchGSEAConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEAConfigSet.md),
[`BenchGeneScoreConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfigSet.md),
[`BenchMemoryConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfigSet.md)

## Examples

``` r
# \donttest{
MethodParamSet(
  Scissor_param = list(alpha = 0.05),
  scPAS_param = list(alpha = 0.05)
)
#> $Scissor_param
#> $Scissor_param$alpha
#> [1] 0.05
#> 
#> 
#> $scAB_param
#> list()
#> 
#> $scPAS_param
#> $scPAS_param$alpha
#> [1] 0.05
#> 
#> 
#> $scPP_param
#> list()
#> 
#> $DEGAS_param
#> list()
#> 
#> $DEGASv2_param
#> list()
#> 
#> $SCIPAC_param
#> list()
#> 
#> $LP_SGL_param
#> list()
#> 
#> $PIPET_param
#> list()
#> 
#> $SIDISH_param
#> list()
#> 
#> $TiRank_param
#> list()
#> 
#> $PhenoMapR_param
#> list()
#> 
# }
```
