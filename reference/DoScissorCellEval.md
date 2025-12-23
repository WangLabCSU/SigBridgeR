# Perform Scissor Cell Evaluation

Evaluates the significance of Scissor-identified cells using benchmark
data and bootstrap resampling.

## Usage

``` r
DoScissorCellEval(
  benchmark_data_path = "path_to_file.RData",
  scissor_res,
  FDR_cutoff = 0.05,
  bootstrap_n = 100L,
  verbose = getFuncOption("verbose"),
  ...
)
```

## Arguments

- benchmark_data_path:

  Path to benchmark data (RData file)

- scissor_res:

  Scissor results from Scissor

- FDR_cutoff:

  FDR threshold for significance (default: `0.05`)

- bootstrap_n:

  Number of bootstrap iterations (default: `100L`)

- verbose:

  Whether to show progress messages (default: `TRUE`)

- ...:

  No usage

## Value

Cell evaluation results with significance assessments

## See also

Other scissor:
[`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md),
[`DoScissorRelTest()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissorRelTest.md)
