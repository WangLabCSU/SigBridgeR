# Screen Function Template

Screen Function Template

## Usage

``` r
DoDEGASv2(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = NULL,
  phenotype_class = c("binary", "survival", "continuous"),
  top_fraction_pos = 0.2,
  sclab = NULL,
  bulk_hvg = TRUE,
  bulk_de = TRUE,
  sc_de = TRUE,
  add_genes = NULL,
  n_hvg = 250L,
  n_bulk_de = 250L,
  n_sc_de = 200L,
  padj.thresh = 0.05,
  only.pos = FALSE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  n_st_classes = length(unique(sc_data$seurat_clusters)),
  loss_type = c("cross_entropy", "log_neg", "rank_loss"),
  transfer_type = c("Wasserstein", "MMD"),
  model_save_dir = "DEGASv2_res",
  lambda1 = 1L,
  lambda2 = 3L,
  lambda3 = 3L,
  tot_seeds = 10L,
  tot_iters = 300L,
  extract_embs = FALSE,
  random_feat = FALSE,
  random_perc = 0.8,
  early_stopping = FALSE,
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

  Phenotype data, either:

  - Named vector (names match `matched_bulk` columns), or

  - Patient survival Data frame with row names matching `matched_bulk`
    columns, colnames named "time" and "status"

- label_type:

  Character specifying phenotype label type

- phenotype_class:

  Type of phenotypic outcome (must be consistent with input data):

  - `"binary"`: Binary traits (e.g., case/control)

  - `"continuous"`: Continuous measurements

  - `"survival"`: Survival infomation

- top_fraction_pos:

  Proportion of cells to be labelled as `"Positive"`

- sclab:

  Metadata for scRNA-seq e.g. cell type, clusters (required), default =
  NULL

- bulk_hvg:

  Whether user wants to use sample (bulk) level high variable genes as
  input in DEGAS. Default is `TRUE`.

- bulk_de:

  Whether user wants to use sample (bulk) differential expression genes
  for phenotype using DEseq2, commonly set Normal or Control as
  reference. Default is `TRUE`.

- sc_de:

  Whether user wants to use single cell differential expression genes
  for each cluster. Default is `TRUE`.

- add_genes:

  The user decides whether to add their own gene list by preference. If
  the user set bulk_hvg, bulk_de and sc_de as FALSE, this gene list is
  necessary. Default is `NULL`.

- n_hvg:

  Number of bulk high variable genes, default = 250.

- n_bulk_de:

  Number of bulk differential expression genes for phenotype, default =
  250.

- n_sc_de:

  Number of genes for each single cell clusters, default = 200.

- padj.thresh:

  Threshold for differential expression analysis, default = 0.05

- only.pos:

  Only return positive markers (FALSE by default)

- min.pct:

  only test genes that are detected in a minimum fraction of min.pct
  cells in either of the two populations. Meant to speed up the function
  by not testing genes that are very infrequently expressed. Default is
  0.01

- logfc.threshold:

  Limit testing to genes which show, on average, at least X-fold
  difference (log-scale) between the two groups of cells. Default is 0.1
  Increasing logfc.threshold speeds up the function, but can miss weaker
  signals. If the `slot` parameter is "scale.data" no filtering is
  performed.

- n_st_classes:

  Integer specifying the number of distinct single-cell/spatial labels
  (cell types, clusters, or spatial annotations). This defines the
  output dimension of the high-resolution (single-cell/spatial) branch,
  and is passed to the Python backend as `high_reso_output_shape`. It is
  required for `Class*` model types, which use single-cell labels to
  supervise the high-resolution branch; for `Blank*` models it is an
  unused dummy value. Typically set to `length(unique(labels))`. Default
  = `NULL`.

- loss_type:

  String specifying the loss function used for training the
  low-resolution (patient) branch. Typically, classification models use
  `cross_entropy`, while survival models can use `log_neg` or
  `rank_loss`. This value must match the expected loss type of the
  chosen `model_type`, otherwise the backend will raise an error.

- transfer_type:

  String indicating the cross-domain alignment method. The backend
  currently supports `Wasserstein` (adversarial with gradient penalty)
  and `MMD` (maximum mean discrepancy). This parameter determines how
  the transfer alignment loss is computed.

- model_save_dir:

  String specifying the output directory for model files. The function
  will automatically create the directory if it does not exist and will
  save configuration files, prediction results, and optionally embedding
  vectors for each submodel.

- lambda1:

  Numeric value, a weighting coefficient that controls the contribution
  of the high-resolution branch (single-cell/spatial supervision loss).
  Only applicable to models with the `Class` prefix. Default = 1.

- lambda2:

  Numeric value, a weighting coefficient that controls the loss strength
  for the low-resolution (patient) branch. This parameter applies to all
  model types. Default = 2.

- lambda3:

  Numeric value, a weighting coefficient that controls the contribution
  of the domain-alignment (transfer) loss. Works jointly with
  `transfer_type` to regulate the alignment between high- and
  low-resolution embedding distributions. Default = 3.

- tot_seeds:

  Integer specifying the number of bagging iterations. Default = 10.

- tot_iters:

  Integer specifying the total number of training iterations (epochs).

- extract_embs:

  Logical value. If `TRUE`, saves the high- and low-resolution hidden
  embeddings as `.npy` files at each save interval for downstream
  visualization or analysis. Default = `FALSE`.

- random_feat:

  Logical value. If `TRUE`, randomly selects a subset of gene features
  before each submodel training. Default = `FALSE`.

- random_perc:

  Numeric value between 0 and 1 specifying the proportion of features to
  retain when `random_feat` is enabled default = 0.8. The same
  proportion is applied to both patient-level and single-cell/spatial
  expression matrices.

- early_stopping:

  Logical value. A reserved flag for early stopping;

- ...:

  Additional arguments passed to the function. Common parameters
  include:

  - verbose: Logical. Whether to print verbose output (default: TRUE).

  - seed: Integer. Random seed for reproducibility (default: 123).

  - Other parameters are passed to
    [`Seurat::FindAllMarkers()`](https://satijalab.org/seurat/reference/FindAllMarkers.html)

## Value

A named list containing:

- scRNA_data:

  Modified single-cell data object with integrated screening results.
