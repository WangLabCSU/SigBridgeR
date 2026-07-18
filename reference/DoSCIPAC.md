# Screen Single-Cell Data Using SCIPAC Algorithm

Performs single-cell phenotype association analysis using the SCIPAC
method. Integrates bulk RNA-seq phenotype information with single-cell
transcriptomic data to identify cell populations associated with
clinical outcomes.

## Usage

``` r
DoSCIPAC(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = "SCIPAC",
  phenotype_class = c("binary", "survival", "continuous"),
  hvg = 1000L,
  do_pca_sc = FALSE,
  n_pc = 60L,
  sc_batch_col = NULL,
  resolution = 2L,
  ela_net_alpha = 0.4,
  bt_size = 50L,
  ncore = 7L,
  ci_alpha = 0.05,
  nfold = 10L,
  ...
)
```

## Arguments

- matched_bulk:

  Matrix or data frame of preprocessed bulk RNA-seq expression data
  (genes × samples). Column names must match names/IDs in `phenotype`.

- sc_data:

  A matrix/Matrix (genes × cells) or a Seurat object containing
  scRNA-seq data to be screened.

- phenotype:

  Phenotype data, either:

  - Named vector (names match `matched_bulk` columns)

  - Data frame with row names matching `matched_bulk` columns,
    containing `"time"` and `"status"` columns for survival analysis

- label_type:

  Character. Context info specifying the screening method identifier.
  Default: `"SCIPAC"`.

- phenotype_class:

  Type of phenotypic outcome:

  `"binary"`

  :   Binary traits (e.g., case/control, responder/non-responder)

  `"continuous"`

  :   Continuous measurements (e.g., age, biomarker levels)

  `"survival"`

  :   Survival information (time-to-event data)

  Partial matching is supported.

- hvg:

  Integer. Number of highly variable genes to use for preprocessing.
  Default: `1000L`.

- do_pca_sc:

  Logical. if `TRUE`, first do PCA on `sc.dat` and use the rotation
  matrix on `bulk.dat`; if `False`, first do PCA on `bulk.data` and use
  the rotation matrix on `sc.dat`. The default is `FALSE`.

- n_pc:

  Integer. Number of principal components to use. Default: `60L`.

- sc_batch_col:

  Character or vector. Batch variable for single-cell data. If
  character, should be a column name in `sc_data` `@metadata`. If
  vector, should be a batch assignment vector matching cell order.
  Default: `NULL` (no batch correction).

- resolution:

  Numeric. Clustering resolution for cell type identification. Higher
  values produce more clusters. Default: `2L`.

- ela_net_alpha:

  Numeric. Elastic net mixing parameter (0 = ridge, 1 = lasso). Default:
  `0.4`.

- bt_size:

  Integer. Bootstrap sample size for stability assessment. Default:
  `50L`.

- ncore:

  Integer. Number of CPU cores for parallel computation. Default: `7L`.

- ci_alpha:

  Numeric. Significance level for confidence intervals. Default: `0.05`.

- nfold:

  Integer. Number of folds for cross-validation for regression models.
  Default: `10L`.

- ...:

  Additional arguments (support `assay`(character), `verbose`(logical) &
  `seed`(integer)).

## Value

A named list containing:

- `scRNA_data`:

  Modified Seurat object with SCIPAC screening results added as metadata
  columns. Includes association statistics, p-values, and confidence
  intervals for each cell cluster.

- `pca_res`:

  PCA rotation results from `SCIPAC::sc.bulk.pca()`

- `cluster_res`:

  Clustering results from `SCIPAC::seurat.ct()`

## Workflow

1.  **Gene overlap**: Identifies common genes between scRNA-seq and bulk
    data

2.  **Preprocessing**: Selects HVGs and preprocesses both datasets

3.  **PCA rotation**: Aligns scRNA-seq and bulk data in shared PCA space

4.  **Clustering**: Identifies cell populations using Seurat clustering

5.  **Association testing**: Tests cluster-phenotype associations using
    elastic net regression with bootstrap validation

6.  **Result integration**: Adds screening results to Seurat object
    metadata

## Requirements

- R package: `SCIPAC` (install from `Exceret/SCIPAC`)

- Gene symbols must be consistent between scRNA-seq and bulk datasets

- Sufficient overlapping genes (typically \>100) for reliable rotation

## Notes

- For survival analysis, `phenotype` must be a data frame with `"time"`
  and `"status"` columns

- Binary phenotypes are automatically converted to factors

- Continuous phenotypes use Gaussian family in elastic net regression

- Screening parameters are stored in `scRNA_data@misc$SCIPAC_params` for
  reproducibility

## See also

`SCIPAC`, `sc.bulk.pca`, `seurat.ct`,
[`RegisterScreenMethod`](https://wanglabcsu.github.io/sigbridger/reference/RegisterScreenMethod.md)

Other screen_method:
[`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md),
[`DoLP_SGL()`](https://wanglabcsu.github.io/sigbridger/reference/DoLP_SGL.md),
[`DoPIPET()`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md),
[`DoSIDISH()`](https://wanglabcsu.github.io/sigbridger/reference/DoSIDISH.md),
[`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md),
[`DoTiRank()`](https://wanglabcsu.github.io/sigbridger/reference/DoTiRank.md),
[`DoscAB()`](https://wanglabcsu.github.io/sigbridger/reference/DoscAB.md),
[`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md),
[`DoscPP()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPP.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Binary phenotype (case/control)
library(SCIPAC)

# Prepare data
bulk_expr <- matrix(rpois(1000, 10), nrow = 200, ncol = 5)
colnames(bulk_expr) <- paste0("Sample", 1:5)
rownames(bulk_expr) <- paste0("Gene", 1:200)

seurat_obj <- CreateSeuratObject(counts = matrix(rpois(2000, 5), nrow = 200, ncol = 10))

phenotype_binary <- c(Sample1 = 1, Sample2 = 0, Sample3 = 1, Sample4 = 0, Sample5 = 1)

# Run SCIPAC screening
result <- DoSCIPAC(
  matched_bulk = bulk_expr,
  sc_data = seurat_obj,
  phenotype = phenotype_binary,
  phenotype_class = "binary",
  hvg = 100,
  resolution = 1,
  verbose = TRUE
)

# Access results
screened_seurat <- result$scRNA_data
head(screened_seurat[[]])

# Example 2: Survival analysis
survival_data <- data.frame(
  time = c(12.5, 24.3, 18.7, 36.2, 30.1),
  status = c(1, 0, 1, 1, 0),
  row.names = paste0("Sample", 1:5)
)

result <- DoSCIPAC(
  matched_bulk = bulk_expr,
  sc_data = seurat_obj,
  phenotype = survival_data,
  phenotype_class = "survival",
  ela_net_alpha = 0.5,
  bt_size = 100,
  ncore = 4
)

# Example 3: Continuous phenotype with custom parameters
continuous_pheno <- c(Sample1 = 25.3, Sample2 = 45.7, Sample3 = 38.2,
                      Sample4 = 52.1, Sample5 = 31.8)

result <- DoSCIPAC(
  matched_bulk = bulk_expr,
  sc_data = seurat_obj,
  phenotype = continuous_pheno,
  phenotype_class = "continuous",
  hvg = 500,
  n_pc = 30,
  resolution = 2,
  ela_net_alpha = 0.3,
  ci_alpha = 0.01,
  seed = 42
)

# View screening parameters used
result$scRNA_data@misc$SCIPAC_params

# Example 4: With batch correction
seurat_obj$batch_id <- sample(c("Batch1", "Batch2"), ncol(seurat_obj), replace = TRUE)

result <- DoSCIPAC(
  matched_bulk = bulk_expr,
  sc_data = seurat_obj,
  phenotype = phenotype_binary,
  phenotype_class = "binary",
  sc_batch_col = "batch_id",
  ncore = 8
)
} # }
```
