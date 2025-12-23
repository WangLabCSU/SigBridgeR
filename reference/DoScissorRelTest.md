# Perform Scissor Reliability Test

Performs reliability testing for Scissor results to assess the stability
and robustness of identified phenotype-associated cells.

## Usage

``` r
DoScissorRelTest(
  scissor_res,
  alpha,
  family,
  cell_num = length(scissor_res$Scissor_pos) + length(scissor_res$Scissor_neg),
  n = 10L,
  nfold = 10L,
  verbose = getFuncOption("verbose"),
  ...
)
```

## Arguments

- scissor_res:

  Scissor results from Scissor

- alpha:

  Alpha parameter used in Scissor, a scalar value between 0 and 1

- family:

  Model family used in Scissor, one of "gaussian", "binomial", "cox"

- cell_num:

  Number of cells identified by Scissor

- n:

  Permutation times (default: `10L`)

- nfold:

  Cross-validation folds (default: `10L`)

- verbose:

  Whether to show progress messages (default: `TRUE`)

- ...:

  No used parameters

## Value

Reliability test results including statistics and p-values

## See also

Other scissor:
[`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md),
[`DoScissorCellEval()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissorCellEval.md)
