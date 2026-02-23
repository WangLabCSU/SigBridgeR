# Run DEGAS Analysis for Single-Cell and Bulk RNA-seq Data Integration

This function performs DEGAS to integrate single-cell and bulk RNA-seq
data, identifying phenotype-associated cells using a bootstrap
aggregated multi-task learning approach.

## Usage

``` r
DoDEGAS(
  select_fraction = 0.05,
  min_thresh = 0.4,
  matched_bulk,
  sc_data,
  phenotype = NULL,
  sc_data.pheno_colname = NULL,
  label_type = "DEGAS",
  phenotype_class = c("binary", "continuous", "survival"),
  tmp_dir = "tmp",
  env_params = list(),
  degas_params = list(),
  normality_test_method = c("jarque-bera", "d'agostino", "kolmogorov-smirnov"),
  ...
)
```

## Arguments

- select_fraction:

  The top percentage of selected cells will be considered as Positive
  cells, without considering how much larger the possible correlation
  coefficient of the observation group is compared to that of the
  control group. Only usedl when `phenotype_class` is "binary" or
  "survival". (default: `0.05`)

- min_thresh:

  DEGAS will calculate the possible correlation coefficients for each
  cell related to the phenotype. When the coefficient of the observation
  group is at least `min_thresh` larger than that of the control group,
  it can be considered related to the phenotype and will be marked as
  Positive. The priority of `min_thresh` is higher than that of
  `select_fraction.` (default: `0.4`)

- matched_bulk:

  Bulk RNA-seq data as matrix or data.frame (rows=genes,
  columns=samples)

- sc_data:

  Single-cell data as Seurat object containing RNA assay

- phenotype:

  Bulk-level phenotype data. For classification: binary matrix with
  one-hot encoding. For survival: matrix with two columns (time and
  event status). Can be NULL, matrix, data.frame, or vector.

- sc_data.pheno_colname:

  Column name for single-cell phenotype in metadata (if available),
  default: `NULL`

- label_type:

  Label type for DEGAS results (default: "`DEGAS"`)

- phenotype_class:

  Type of phenotype: "binary" (classification), "continuous", or
  "survival"

- tmp_dir:

  Temporary directory for intermediate files (default: `"tmp"`)

- env_params:

  List of environment parameters for Python setup including:

  - env.name: environment name (default: `"r-reticulate-degas"`)

  - env.type: environment type "conda", "environment", or "venv"
    (default: `"environment"`)

  - env.method: environment setup method "system", "conda" (default:
    `"system"`)

  - env.file: path to environment file (default:
    system.file("conda/DEGAS_environment.yml", package = "SigBridgeR"))

  - env.python_version: Python version (default: "3.9.15")

  - env.packages: named vector of Python packages and versions (default:
    c("tensorflow" = "2.4.1", "protobuf" = "3.20" ,"numpy" = "any"))

  - env.recreate: whether to recreate environment (default: `FALSE`)

  - env.use_conda_forge: whether to use conda-forge channel (conda only,
    default: `TRUE`)

  - env.verbose: verbose output (default: `FALSE`)

- degas_params:

  List of DEGAS algorithm parameters including:

  - DEGAS.model_type: model type ("BlankClass", "ClassBlank",
    "ClassClass", "ClassCox", "BlankCox")

  - DEGAS.architecture: "Standard" (feed forward) or "DenseNet" (dense
    net), default: `"DenseNet"`

  - DEGAS.ff_depth: number of layers in model (\>=1, default: `3`)

  - DEGAS.pyloc: path to Python executable (default: `NULL`, automatic
    detection)

  - DEGAS.bag_depth: bootstrap aggregation depth (\>=1, default: `5`)

  - DEGAS.train_steps: training steps (default: `2000`)

  - DEGAS.scbatch_sz: single-cell batch size (default: `200`)

  - DEGAS.patbatch_sz: patient batch size (default: `50`)

  - DEGAS.hidden_feats: hidden features (default: `50`)

  - DEGAS.do_prc: dropout percentage (default: `0.5`)

  - DEGAS.lambda1: regularization parameter 1 (default: `3.0`)

  - DEGAS.lambda2: regularization parameter 2 (default: `3.0`)

  - DEGAS.lambda3: regularization parameter 3 (default: `3.0`)

  - DEGAS.seed: random seed (default: `2`)

- normality_test_method:

  Method for normality testing: "jarque-bera", "d'agostino", or
  "kolmogorov-smirnov"

- ...:

  Additional arguments. Currently supports:

  - `verbose`: Logical indicating whether to print progress messages.
    Defaults to `TRUE`.

  - `assay`: Name of assay to use. Defaults to "RNA".

## Value

A list containing:

- scRNA_data: Seurat object with DEGAS labels added to metadata

- model: The model trained using the input data, andit can be used for
  cell classification prediction.

- DEGAS_prediction: Data table with DEGAS predictions containing:

  - Predicted label probabilities for each cell

  - Cell labels ("Positive"/"Other") based on selection criteria

  - Difference scores for binary phenotypes

  - Cell identifiers

## Details

The function performs the following steps:

1.  Validates input data and parameters

2.  Sets up Python environment with required dependencies

3.  Trains bootstrap aggregated DEGAS model using `runCCMTLBag`

4.  Generates cell-level predictions using `predClassBag`

5.  Applies statistical testing to identify phenotype-associated cells

6.  Labels cells as "Positive" or "Other" based on selection criteria

Model type is automatically determined:

- BlankClass: only bulk phenotype specified (scLab = NULL)

- ClassBlank: only single-cell phenotype specified (patLab = NULL)

- ClassClass: both single-cell and bulk phenotypes specified

- ClassCox: single-cell phenotype + bulk survival data

- BlankCox: only bulk survival data specified

## References

Johnson TS, Yu CY, Huang Z, Xu S, Wang T, Dong C, et al. Diagnostic
Evidence GAuge of Single cells (DEGAS): a flexible deep transfer
learning framework for prioritizing cells in relation to disease. Genome
Med. 2022 Feb 1;14(1):11.

## See also

Other screen_method:
[`DoLP_SGL()`](https://wanglabcsu.github.io/sigbridger/reference/DoLP_SGL.md),
[`DoPIPET()`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md),
[`DoSIDISH()`](https://wanglabcsu.github.io/sigbridger/reference/DoSIDISH.md),
[`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md),
[`DoscAB()`](https://wanglabcsu.github.io/sigbridger/reference/DoscAB.md),
[`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md),
[`DoscPP()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPP.md)

Other DEGAS:
[`DEGASEnvSet()`](https://wanglabcsu.github.io/sigbridger/reference/DEGASEnvSet.md),
[`DEGASModelDetect()`](https://wanglabcsu.github.io/sigbridger/reference/DEGASModelDetect.md),
[`DEGASParamSet()`](https://wanglabcsu.github.io/sigbridger/reference/DEGASParamSet.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Binary classification example
result <- DoDEGAS(
  select_fraction = 0.05, # `select_fraction` only used in binary and survival phenotyping
  matched_bulk = bulk_matrix,
  sc_data = seurat_obj,
  phenotype = bulk_phenotype,
  phenotype_class = "binary"
)

# Survival analysis example
result <- DoDEGAS(
  select_fraction = 0.05, # `select_fraction` only used in binary and survival phenotyping
  matched_bulk = bulk_matrix,
  sc_data = seurat_obj,
  phenotype = survival_data,
  phenotype_class = "survival"
)
} # }
```
