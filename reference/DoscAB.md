# Perform scAB Screening Analysis

Implements the scAB algorithm to identify phenotype-associated cell
subpopulations in single-cell RNA-seq data by integrating matched bulk
expression and phenotype information. Uses non-negative matrix
factorization (NMF) with dual regularization for phenotype association
and cell-cell similarity.

## Usage

``` r
DoscAB(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = "scAB",
  phenotype_class = c("binary", "survival"),
  alpha = c(0.005, NULL),
  alpha_2 = c(0.005, NULL),
  maxiter = 2000L,
  tred = 2L,
  ...
)
```

## Arguments

- matched_bulk:

  Normalized bulk expression matrix (genes × samples) where: - Columns
  match `phenotype` row names - Genes match features in `sc_data`

- sc_data:

  Seurat object containing preprocessed single-cell data:

- phenotype:

  Data frame with clinical annotations where: - Rows correspond to
  `matched_bulk` columns - For survival: contains `time` and `status`
  columns

- label_type:

  Character specifying phenotype label type (e.g., "SBS1", "time"),
  stored in `scRNA_data@misc`

- phenotype_class:

  Analysis mode: - `"binary"`: Case-control design (e.g.,
  responder/non-responder) - `"survival"`: Time-to-event analysis
  data.frame

- alpha:

  Coefficient of phenotype regularization (default=`0.005`).

- alpha_2:

  Coefficent of cell-cell similarity regularization (default=`0.005`).

- maxiter:

  Maximum number of iterations for NMF (default=2000).

- tred:

  Z-score threshold in finding subsets (default=2).

- ...:

  Additional arguments. Currently supports:

  - `verbose`: Logical indicating whether to print progress messages.
    Defaults to `TRUE`.

  - `seed`: For reproducibility, default is `123L`

  - `parallel`: Logical indicating whether to use parallel processing.
    Defaults to `FALSE`.

  - Other arguments are passed to
    [`scAB::create_scAB.v5()`](https://rdrr.io/pkg/scAB/man/create_scAB.v5.html)

  - `assay`: Character specifying the assay to use. Defaults to `"RNA"`.

## Value

A list containing:

- scRNA_data:

  Filtered Seurat object with selected cells

- scAB_result:

  scAB screening result

## LICENSE

Licensed under the GNU General Public License version 3 (GPL-3.0). A
copy of the license is available at
<https://www.gnu.org/licenses/gpl-3.0.en.html>.

## References

Zhang Q, Jin S, Zou X. scAB detects multiresolution cell states with
clinical significance by integrating single-cell genomics and bulk
sequencing data. Nucleic Acids Research. 2022 Nov 28;50(21):12112–30.

## See also

Other screen_method:
[`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md),
[`DoLP_SGL()`](https://wanglabcsu.github.io/sigbridger/reference/DoLP_SGL.md),
[`DoPIPET()`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md),
[`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md),
[`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md),
[`DoscPP()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPP.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Binary phenotype example
result <- DoscAB(
  matched_bulk = bulk_matrix,
  sc_data = seurat_obj,
  phenotype = clinical_df,
  label_type = "disease_status",
  phenotype_class = "binary",
  alpha = 0.005,
  alpha_2 = 0.005,
  maxiter = 2000,
  tred = 2
)
} # }
```
