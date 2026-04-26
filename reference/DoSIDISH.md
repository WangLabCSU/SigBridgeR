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
  sidish_params = list(),
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

- sidish_params:

  List of SIDISH algorithm parameters including: **Preprocessing
  parameters:**

  - `patient_id`: column name for patient identifier in metadata
    (default: `"Sample"`)

  - `celltype_name`: column name for cell type annotation in metadata
    (default: `"celltype_major"`)

  - `processed`: whether input data is already preprocessed (default:
    `TRUE`)

  - `n_genes_by_counts`: minimum number of genes expressed per cell
    filter threshold (default: `5000`)

  - `pct_counts_mt`: maximum percentage of mitochondrial genes filter
    threshold (default: `10`)

  - `batch_correction`: whether to perform batch correction (default:
    `FALSE`)

  - `survival_`: column name for survival time in phenotype data
    (default: `"time"`)

  - `status`: column name for event status in phenotype data (default:
    `"status"`)

  **Execution environment:**

  - `device`: computation device, `"cuda"` for GPU acceleration or
    `"cpu"` for CPU-only (default: `"cuda"`)

  - `use_spatial_graph`: whether to use spatial graph information
    (default: `FALSE`)

  - `k_neighbors`: number of neighbors for graph construction (default:
    `NULL`, auto-detected)

  **Phase 1: VAE training parameters:**

  - `phase1_epochs`: total epochs for VAE training (default: `225`)

  - `phase1_i_epochs`: interval epochs for VAE intermediate evaluation
    (default: `20`)

  - `phase1_latent_size`: dimensionality of latent space (default: `32`)

  - `phase1_layer_dims`: hidden layer dimensions as integer vector
    (default: `c(512, 128)`)

  - `phase1_batch_size`: batch size for VAE training (default: `256`)

  - `phase1_optimizer`: optimizer algorithm (default: `"Adam"`)

  - `phase1_lr`: learning rate for VAE encoder/decoder (default: `1e-4`)

  - `phase1_lr_3`: learning rate for additional VAE component (default:
    `1e-4`)

  - `phase1_dropout`: dropout rate for VAE layers (default: `0`)

  - `phase1_type`: VAE layer type, `"Dense"` or `"Normal"` (default:
    `"Dense"`)

  **Phase 2: Deep Cox training parameters:**

  - `phase2_epochs`: total epochs for Cox model training (default:
    `500`)

  - `phase2_hidden`: number of hidden units in Cox model (default:
    `128`)

  - `phase2_lr`: learning rate for Cox model (default: `1e-4`)

  - `phase2_dropout`: dropout rate for Cox model (default: `0`)

  - `phase2_test_size`: proportion of data held out for testing
    (default: `0.2`)

  - `phase2_batch_size_bulk`: batch size for bulk data in Cox training
    (default: `256`)

  **Training & risk definition parameters:**

  - `train_iterations`: number of risk score iteration rounds (default:
    `5`)

  - `train_percentile`: percentile threshold for high-risk cell
    selection (default: `0.95`)

  - `train_steepness`: steepness parameter for risk score transformation
    (default: `30`)

  - `train_path`: directory path for saving intermediate results
    (default: `"./SIDISH_res/"`)

  - `train_num_workers`: number of data loading workers (default: `0`)

  - `train_distribution_fit`: distribution fitting method, `"fitted"` or
    `"default"` (default: `"fitted"`)

- env_params:

  List of environment parameters for Python setup including:

  - `env.name`: conda/environment name (default:
    `"r-reticulate-sidish-nvidia"` for CUDA or
    `"r-reticulate-sidish-cpu"` for CPU)

  - `env.type`: environment type, `"conda"`, `"environment"`, or
    `"venv"` (default: `"conda"`)

  - `env.method`: environment setup method, `"system"` or `"conda"`
    (default: `"environment"`)

  - `env.file`: path to environment YAML file (default:
    `system.file("conda/SIDISH_nvidia_environment.yml", package = "SigBridgeR")`
    or CPU variant)

  - `env.python_version`: Python version (default: `"3.12.12"`)

  - `env.packages`: named vector of Python packages and versions
    (default: `c("numpy" = "1.26.4")`, more packages included via
    env.file)

  - `env.recreate`: whether to recreate the environment if it already
    exists (default: `FALSE`)

  - `env.use_conda_forge`: whether to use the conda-forge channel (conda
    only, default: `TRUE`)

  - `env.verbose`: verbose output during environment setup (default:
    value from `getFuncOption("verbose")`)

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
