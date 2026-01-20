# Perform scPAS Screening Analysis

This function performs scPAS screening analysis by integrating bulk and
single-cell RNA-seq data. It includes data filtering steps and wraps the
core scPAS::scPAS function.

## Usage

``` r
DoscPAS(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = "scPAS",
  assay = "RNA",
  imputation = FALSE,
  imputation_method = c("KNN", "ALRA"),
  nfeature = 3000L,
  alpha = c(0.01, NULL),
  cutoff = 0.2,
  network_class = c("SC", "bulk"),
  family = c("cox", "gaussian", "binomial"),
  permutation_times = 2000L,
  FDR_threshold = 0.05,
  independent = TRUE,
  ...
)
```

## Arguments

- matched_bulk:

  Bulk RNA-seq data (genes x samples)

- sc_data:

  Single-cell RNA-seq data (Seurat object and preprocessed)

- phenotype:

  Phenotype data frame with sample annotations

- label_type:

  Character specifying phenotype label type (e.g., "SBS1", "time"),
  stored in `scRNA_data@misc`. (default: "scPAS")

- assay:

  Assay to use from sc_data (default: 'RNA')

- imputation:

  Logical, whether to perform imputation (default: FALSE)

- imputation_method:

  Character. Name of alternative method for imputation. (options: "KNN",
  "ALRA")

- nfeature:

  Number of features to select (default: 3000, indicating that the top
  3000 highly variable genes are selected for model training

- alpha:

  Numeric. Significance threshold. Parameter used to balance the effect
  of the l1 norm and the network-based penalties. It can be a number or
  a searching vector. If alpha = NULL, a default searching vector is
  used. The range of alpha is in `[0,1]`. A larger alpha lays more
  emphasis on the l1 norm. (default: 0.01)

- cutoff:

  Numeric. Cutoff value for selecting the optimal alpha value when alpha
  = NULL. (default: 0.2)

- network_class:

  Network class to use (default: 'SC', indicating gene-gene similarity
  networks derived from single-cell data. The other one is 'bulk'.)

- family:

  Model family for analysis (options: "cox", "gaussian", "binomial")

- permutation_times:

  Number of permutations to perform (default: 2000)

- FDR_threshold:

  Numeric. FDR value threshold for identifying phenotype-associated
  cells (default: 0.05)

- independent:

  Logical. The background distribution of risk scores is constructed
  independently of each cell. (default: TRUE)

- ...:

  Additional arguments. Currently supports:

  - `verbose`: Logical indicating whether to print progress messages.
    Defaults to `TRUE`.

  - `seed`: For reproducibility, default is `123L`

## Value

A Seurat object from scPAS analysis

## LICENSE

Licensed under the GNU General Public License version 3 (GPL-3.0). A
copy of the license is available at
<https://www.gnu.org/licenses/gpl-3.0.en.html>.

## References

Xie A, Wang H, Zhao J, Wang Z, Xu J, Xu Y. scPAS: single-cell
phenotype-associated subpopulation identifier. Briefings in
Bioinformatics. 2024 Nov 22;26(1):bbae655.

## See also

Other screen_method:
[`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md),
[`DoLP_SGL()`](https://wanglabcsu.github.io/sigbridger/reference/DoLP_SGL.md),
[`DoPIPET()`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md),
[`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md),
[`DoscAB()`](https://wanglabcsu.github.io/sigbridger/reference/DoscAB.md),
[`DoscPP()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPP.md)
