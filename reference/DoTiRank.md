# Perform TiRank Screening Analysis

Integrates single-cell and bulk RNA-seq data using the TiRank deep
learning framework to identify phenotype-associated cells. TiRank
employs a rank-based approach with neural network models to score and
classify cells based on their association with the phenotype of
interest.

## Usage

``` r
DoTiRank(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = "TiRank",
  phenotype_class = c("binary", "survival", "continuous"),
  tirank_params = list(validation_proportion = 0.15, sampling_thresh = 0.5, sampling_mode
    = c("smote", "downsample", "upsample", "tomeklinks"), top_var_genes = 2000L,
    top_gene_pairs = 1000L, p_value_threshold = 0.05, max_cutoff = 0.8, min_cutoff =
    -0.8, nhead = 2L, nhid1 = 96L, nhid2 = 8L, n_output = 32L, nlayers = 3L, n_pred = 2L,
    dropout = 0.5, encoder_type = c("MLP", "Transformer", "DenseNet"), infer_mode =
    c("SC", "ST"), n_trials = 5L, do_reject = TRUE, tolerance = 0.05),
  save_path = "./TiRank_res",
  load_cache = NULL,
  ...
)
```

## Arguments

- matched_bulk:

  Matrix or data frame of preprocessed bulk RNA-seq expression data
  (genes x samples). Column names must match names/IDs in `phenotype`.

- sc_data:

  A matrix/Matrix (genes x cells) or a Seurat object containing
  scRNA-seq data to be screened.

- phenotype:

  Phenotype data, either: - Named vector (names match `matched_bulk`
  columns), or - Patient survival data frame with row names matching
  `matched_bulk` columns, colnames named "time" and "status"

- label_type:

  Character specifying phenotype label type (default: `"TiRank"`).

- phenotype_class:

  Type of phenotypic outcome (must be consistent with input data): -
  `"binary"`: Binary traits (e.g., case/control) - `"continuous"`:
  Continuous measurements - `"survival"`: Survival information

- tirank_params:

  List of TiRank algorithm parameters:

  validation_proportion

  :   Proportion of bulk data held out for validation (default: `0.15`).

  sampling_thresh

  :   Threshold for resampling strategy (default: `0.5`).

  sampling_mode

  :   Resampling method: `"smote"`, `"downsample"`, `"upsample"`, or
      `"tomeklinks"`.

  top_var_genes

  :   Number of top variable genes to select (default: `2000L`).

  top_gene_pairs

  :   Number of top gene pairs for ranking (default: `1000L`).

  p_value_threshold

  :   P-value threshold for feature selection (default: `0.05`).

  max_cutoff

  :   Upper cutoff for correlation filtering (default: `0.8`).

  min_cutoff

  :   Lower cutoff for correlation filtering (default: `-0.8`).

  nhead

  :   Number of attention heads (default: `2L`).

  nhid1

  :   Size of first hidden layer (default: `96L`).

  nhid2

  :   Size of second hidden layer (default: `8L`).

  n_output

  :   Output dimension of encoder (default: `32L`).

  nlayers

  :   Number of encoder layers (default: `3L`).

  n_pred

  :   Number of prediction heads (default: `2L`).

  dropout

  :   Dropout rate for regularization (default: `0.5`).

  encoder_type

  :   Encoder architecture: `"MLP"`, `"Transformer"`, or `"DenseNet"`.

  infer_mode

  :   Inference mode: `"SC"` (single-cell) or `"ST"` (spatial
      transcriptomics).

  n_trials

  :   Number of repeated training trials (default: `5L`).

  do_reject

  :   Whether to apply rejection criteria to uncertain predictions
      (default: `TRUE`).

  tolerance

  :   Tolerance threshold for rejection (default: `0.05`).

- save_path:

  (Soft-deprecated) Directory path for saving intermediate and final
  results (default: `"./TiRank_res"`). Acts as fallback for `load_cache`
  and `save_cache` when not specified via `...`. Prefer using
  `load_cache` / `save_cache` in `...` instead. See
  [`CacheSetHere()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSetHere.md).

- load_cache:

  (Soft-deprecated) Optional path to cached data (default: `NULL`).
  Prefer using `load_cache` in `...` instead.

- ...:

  Additional arguments passed to the function. Common parameters
  include:

  verbose

  :   Logical. Whether to print verbose output (default: TRUE).

  seed

  :   Integer. Random seed for reproducibility.

  assay

  :   Character. Name of assay to use from Seurat object (default:
      `"RNA"`).

  load_cache

  :   Cache directory path for loading cached data. Supports root-level,
      cache-level, or parent-level paths. See
      [`CacheSetHere()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSetHere.md).

  save_cache

  :   Cache directory path for saving results. Supports root-level or
      parent-level paths. See
      [`CacheSetHere()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSetHere.md).

## Value

A named list containing:

- scRNA_data:

  Modified single-cell data object with integrated screening results
  added as metadata, including `TiRank_Reject`, `TiRank_Rank_Score`, and
  `TiRank` (Positive/Neutral/Negative) columns, plus `TiRank_para` and
  `TiRank_type` stored in misc slot.

- cell_cell_distance:

  Computed cell-cell similarity/distance matrix.

## Details

The TiRank screening workflow consists of the following steps:

1.  **Data preprocessing:** Normalizes bulk expression data and
    validates compatibility with phenotype information.

2.  **Expression transfer:** Transfers the single-cell expression
    profile into the TiRank-compatible format.

3.  **Validation set generation:** Splits bulk data into training and
    validation sets according to `validation_proportion`.

4.  **Resampling:** Applies the specified sampling strategy to address
    class imbalance in bulk data.

5.  **Cell-cell similarity:** Computes cell-cell distance or similarity
    matrix from single-cell data.

6.  **Model training:** Trains a TiRank neural network (MLP,
    Transformer, or DenseNet encoder) using bulk expression and
    phenotype.

7.  **Screening:** Applies the trained model to score each cell,
    producing rank-based predictions with rejection handling.

8.  **Label assignment:** Classifies cells as `"Positive"`, `"Neutral"`,
    or `"Negative"` based on rank scores.

## See also

Other screen_method:
[`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md),
[`DoLP_SGL()`](https://wanglabcsu.github.io/sigbridger/reference/DoLP_SGL.md),
[`DoPIPET()`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md),
[`DoSCIPAC()`](https://wanglabcsu.github.io/sigbridger/reference/DoSCIPAC.md),
[`DoSIDISH()`](https://wanglabcsu.github.io/sigbridger/reference/DoSIDISH.md),
[`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md),
[`DoscAB()`](https://wanglabcsu.github.io/sigbridger/reference/DoscAB.md),
[`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md),
[`DoscPP()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPP.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Binary classification example
result <- DoTiRank(
  matched_bulk = bulk_matrix,
  sc_data = seurat_obj,
  phenotype = bulk_phenotype,
  phenotype_class = "binary",
  save_path = "./TiRank_res"
)

# Survival analysis example
result <- DoTiRank(
  matched_bulk = bulk_matrix,
  sc_data = seurat_obj,
  phenotype = survival_data,
  phenotype_class = "survival",
  tirank_params = list(
    encoder_type = "Transformer",
    n_trials = 3L
  ),
  save_path = "./TiRank_survival_res"
)

# Access results
modified_seurat <- result$scRNA_data
head(modified_seurat[[]])
} # }
```
