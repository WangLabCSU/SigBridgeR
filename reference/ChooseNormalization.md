# Data-Driven Selection of Single-Cell Normalization Methods

A quantitative framework for selecting optimal normalization strategies
(e.g., SCTransform vs. LogNormalization) based on diagnostic metrics
rather than heuristics. Evaluates three critical dimensions of
preprocessing quality:

1.  **Variance stabilization**: Decoupling of mean-variance relationship
    in normalized expression (lower correlation = better).

2.  **Biological signal retention**: Preservation of known marker genes
    within highly variable genes (higher retention = better).

3.  **Dropout robustness**: Removal of technical dropout bias from
    normalized values (lower correlation with dropout rate = better).

Methods are ranked using a weighted composite score. Designed for
head-to-head comparison of preprocessed Seurat objects.

## Usage

``` r
ChooseNormalization(
  ...,
  subset_size = integer(),
  known_hvgs = list(),
  n_hvgs = 2000L,
  low_expressed_thresh = 0.2,
  weight = c(variance_stability = 0.4, marker_signal = 0.35, dropout_robustness = 0.25)
)
```

## Arguments

- ...:

  Named arguments where each value is a `Seurat` object representing a
  distinct preprocessing strategy (e.g.,
  `SCT = sct_obj, Log = log_obj`). **Requirements**:

  - Must be named (names become method identifiers)

  - Must contain normalized data in the `data` slot of the specified
    assay/layer

  - Must have identical cell counts (for fair comparison)

- subset_size:

  Integer. Number of cells to subsample for diagnostics. If missing or
  empty, defaults to `min(10000, total_cells)`. Smaller subsets
  accelerate computation with minimal accuracy loss for large datasets
  (\>50k cells).

- known_hvgs:

  Named list of canonical marker genes for cell types. Format:
  `list(cell_type1 = c("geneA", "geneB"), cell_type2 = c("geneC", ...))`.
  Used to compute `marker_retention` metric. If `NULL` (default), this
  metric is omitted from scoring.

- n_hvgs:

  Integer. Number of top highly variable genes to evaluate for marker
  retention. Default: `2000L`. Only relevant when `known_hvgs` is
  provided.

- low_expressed_thresh:

  Numeric (0–1). Quantile threshold for filtering lowly expressed genes
  during variance-mean correlation calculation. Genes below this
  quantile of mean expression are excluded. Default: `0.2` (bottom 20%
  excluded).

- weight:

  Named numeric vector specifying weights for composite scoring. Must
  contain exactly these components summing to 1:

  `variance_stability`

  :   Weight for variance-mean decoupling (default: 0.4)

  `marker_signal`

  :   Weight for marker gene retention (default: 0.35)

  `dropout_robustness`

  :   Weight for dropout bias removal (default: 0.25)

## Value

A list containing:

- `metrics`:

  A `data.table` with diagnostic metrics per method:

  - `variance_mean_cor`: Pearson correlation between log10(mean) and
    log10(variance) of normalized expression (lower = better)

  - `marker_retention`: Proportion of known markers in top HVGs (higher
    = better; NA if `known_hvgs` not provided)

  - `mean_dropout_residual`: Absolute Spearman correlation between
    dropout rate (from counts) and normalized means (lower = better)

  - `composite_score`: Weighted combination of normalized metrics (0–1
    scale)

  - `rank`: Method ranking (1 = best)

- `recommendation`:

  Character string naming the top-ranked method

- `plots`:

  (Optional) Diagnostic visualizations if `ggplot2` available

## Workflow

1.  Validates input objects (naming, cell count consistency, data slot
    presence)

2.  Subsamples cells (if needed) for computational efficiency

3.  Computes three core metrics per method:

    - Variance stabilization:

      Correlation between log-transformed mean and variance of
      normalized expression

    - Marker retention:

      Overlap between user-provided markers and top HVGs

    - Dropout robustness:

      Correlation between gene dropout rate (from counts) and normalized
      expression means

4.  Normalizes metrics to 0-1 scale (inverting where lower=better)

5.  Computes weighted composite score and ranks methods

6.  Returns top recommendation with full diagnostic report

## Notes

- Requires *pre-normalized* Seurat objects—this function **does not**
  perform normalization itself. Users must first run `SCTransform()`,
  `NormalizeData()`, etc., and store results in the `data` slot.

- Marker retention metric is only computed when `known_hvgs` is
  provided. Without it, scoring relies solely on variance stabilization
  and dropout robustness.

- Weights should sum to 1 (not enforced but recommended for
  interpretable scores).

## Examples

``` r
if (FALSE) { # \dontrun{
# Compare two preprocessed Seurat objects
sct_obj <- SCTransform(seurat_raw, verbose = FALSE)
log_obj <- NormalizeData(seurat_raw, normalization.method = "LogNormalize")

result <- ChooseNormalization(
  SCT = sct_obj,
  Log = log_obj,
  subset_size = 5000,
  verbose = TRUE
)

# Top recommendation
result$recommendation

# Full metrics table
result$metrics[, .(method, composite_score, rank)]

# With marker genes for refined scoring
markers <- list(
  T_cell = c("CD3D", "CD3E", "CD8A"),
  B_cell = c("CD79A", "MS4A1", "CD19"),
  Myeloid = c("CD14", "LYZ", "FCGR3A")
)

result <- ChooseNormalization(
  SCT = sct_obj,
  Log = log_obj,
  known_hvgs = markers,
  weight = c(variance_stability = 0.3, marker_signal = 0.5, dropout_robustness = 0.2)
)
} # }
```
