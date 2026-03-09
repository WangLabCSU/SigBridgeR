# Perform SIDISH Screening Analysis

Perform SIDISH Screening Analysis

## Usage

``` r
DoSIDISH(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = NULL,
  phenotype_class = "survival",
  assay = "RNA",
  sidish_param = list(),
  env_params = list(),
  ...
)
```

## Arguments

- matched_bulk:

  Matrix or data frame of preprocessed bulk RNA-seq expression data
  (genes x samples). Column names must match names/IDs in `phenotype`.

- sc_data:

  A Seurat object containing scRNA-seq data to be screened.

- phenotype:

  Phenotype data, either: - Patient survival Data frame with row names
  matching `matched_bulk` columns, colnames named "time" and "status"

- label_type:

  Character specifying phenotype label type

- phenotype_class:

  Type of phenotypic outcome (must be consistent with input data): -
  `"survival"`: Survival infomation

- assay:

  Seurat assay name, default: `"RNA"`.

- sidish_param:

  Parameters adjusting SIDISH. Use `SigBridgeR:::SIDISHParamSet()` to
  see all available parameters and fetch default values.

- env_params:

  Parameters adjusting python environment. Use
  `SigBridgeR:::SIDISHEnvSet()`to see all available parameters and fetch
  default values.

- ...:

  Additional arguments passed to the function. Common parameters
  include:

  verbose

  :   Logical. Whether to print verbose output (default: `TRUE`).

## Value

A named list containing:

- scRNA_data:

  Modified single-cell data object with integrated screening results.

## See also

Other screen_method:
[`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md),
[`DoLP_SGL()`](https://wanglabcsu.github.io/sigbridger/reference/DoLP_SGL.md),
[`DoPIPET()`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md),
[`DoSCIPAC()`](https://wanglabcsu.github.io/sigbridger/reference/DoSCIPAC.md),
[`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md),
[`DoscAB()`](https://wanglabcsu.github.io/sigbridger/reference/DoscAB.md),
[`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md),
[`DoscPP()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPP.md)
