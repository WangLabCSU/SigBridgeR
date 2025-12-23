# Perform LP-SGL Screening Analysis

Identifies phenotype-associated cell subpopulations using
Lasso-Penalized Sparse Group Lasso (LP-SGL) with Leiden community
detection. This method integrates bulk and single-cell RNA-seq data to
identify cell subpopulations associated with phenotypic outcomes.

## Usage

``` r
DoLP_SGL(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = "LP_SGL",
  LPSGL_family = c("logit", "cox", "linear"),
  resolution = 0.6,
  alpha = 0.5,
  nfold = 5,
  dge_analysis = list(run = FALSE, logFC_threshold = 1, pval_threshold = 0.05),
  ...
)
```

## Arguments

- matched_bulk:

  Bulk expression matrix (features × samples)

- sc_data:

  Single-cell RNA-seq data (Seurat object)

- phenotype:

  Binary phenotype vector for bulk samples

- label_type:

  Character specifying phenotype label type (default: "LP_SGL")

- LPSGL_family:

  Type of regression model: "`logit`" (logistic), "`cox`" (Cox), or
  "`linear`" (linear regression)

- resolution:

  Resolution parameter for Leiden clustering (default: `0.6`)

- alpha:

  Alpha parameter for SGL balancing L1 and L2 penalties (default: `0.5`)

- nfold:

  Number of folds for cross-validation (default: `5`)

- dge_analysis:

  List controlling differential expression analysis:

  - run: Whether to run DEG analysis (default: `FALSE`)

  - logFC_threshold: Log fold change threshold (default: `1`)

  - pval_threshold: P-value threshold (default: `0.05`)

- ...:

  Additional arguments passed to preprocessing functions, e.g.:

  - verbose: Whether to print progress messages (default: `TRUE`)

  - seed: Random seed for reproducibility (default: `123L`)

## Value

A list containing:

- scRNA_data:

  Seurat object with LP-SGL results integrated

- sgl_fit:

  Fitted SGL model object

- cvfit:

  Cross-validation results

- dge_res:

  Differential expression results if requested (NULL otherwise)

## References

Li J, Zhang H, Mu B, Zuo H, Zhou K. Identifying phenotype-associated
subpopulations through LP_SGL. Briefings in Bioinformatics. 2023 Nov
22;25(1):bbad424.

## See also

Other screen_method:
[`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md),
[`DoPIPET()`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md),
[`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md),
[`DoscAB()`](https://wanglabcsu.github.io/sigbridger/reference/DoscAB.md),
[`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md),
[`DoscPP()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPP.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Example using simulated data
set.seed(123)

# Create simulated data
bulk_data <- matrix(rnorm(1000*50), nrow=1000, ncol=50)
sc_data <- matrix(rnorm(1000*500), nrow=1000, ncol=500)
phenotype <- rep(c(0, 1), each=25)

# Run LP-SGL analysis
results <- DoLP_SGL(
matched_bulk = bulk_data,
sc_data = sc_data,
phenotype = phenotype,
LPSGL_family = "logit",
resolution = 0.6,
dge_analysis = list(run = TRUE, logFC_threshold = 1, pval_threshold = 0.05)
)

# Access results
lpsgl_seurat <- results$scRNA_data
sgl_model <- results$sgl_fit
deg_results <- results$dge_res
} # }
```
