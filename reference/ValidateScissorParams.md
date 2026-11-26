# Validate and Prepare DoScissor Parameters

Internal helper that handles deprecated argument warnings, input
validation, dots processing, and default-value resolution for
[`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md).

## Usage

``` r
ValidateScissorParams(
  phenotype,
  label_type,
  phenotype_class = c("binary", "continuous", "survival"),
  alpha,
  cutoff,
  family = lifecycle::deprecated(),
  reliability_test,
  cell_evaluation,
  path2load_scissor_cache = lifecycle::deprecated(),
  path2save_scissor_inputs = lifecycle::deprecated(),
  ...
)
```

## Arguments

- phenotype, label_type, alpha, cutoff, family:

  Forwarded from
  [`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md).

- reliability_test, cell_evaluation:

  Forwarded from
  [`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md).

- path2load_scissor_cache, path2save_scissor_inputs:

  Deprecated arguments.

- ...:

  Additional dots forwarded from
  [`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md).

## Value

A named list with elements: `family`, `verbose`, `seed`, `assay`,
`load_cache`, `save_cache`, `reliability_test`, `cell_evaluation`,
`label_type`, `phenotype_class_cache`, `dots`.

## See also

Other scissor:
[`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md),
[`DoScissorCellEval()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissorCellEval.md),
[`DoScissorRelTest()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissorRelTest.md)
