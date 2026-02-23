# Perform Scissor Screening Analysis

Identifies phenotype-associated cell subpopulations in single-cell data
using regularized regression on matched bulk expression profiles.
Scissor integrates bulk and single-cell RNA-seq data to identify cells
that are significantly associated with phenotypic outcomes.

## Usage

``` r
DoScissor(
   path2load_scissor_cache = NULL,
   path2save_scissor_inputs = "Scissor_inputs.RData",
   matched_bulk,
   sc_data,
   phenotype,
   label_type = "scissor",
   alpha = c(0.05, NULL),
   cutoff = 0.2,
   family = c("gaussian", "binomial", "cox"),
   reliability_test = list(
     run = FALSE, # whether to run reliability test
     n = 10L, # permutation times
     nfold = 10L # cross validation folds
   ),
   cell_evaluation = list(
     run = FALSE, # whether to run cell evaluation
     benchmark_data = "path_to_file.RData", # path to benchmark data
     FDR_cutoff = 0.05,
     bootstrap_n = 100L
   ),
   ...
)
```

## Arguments

- path2load_scissor_cache:

  Path to precomputed Scissor inputs (RData file). If provided, skips
  recomputation (default: NULL).

- path2save_scissor_inputs:

  Path to save intermediate files (default: "Scissor_inputs.RData").

- matched_bulk:

  Normalized bulk expression matrix (features × samples). Column names
  must match `phenotype` identifiers.

- sc_data:

  Seurat object containing single-cell RNA-seq data.

- phenotype:

  Clinical outcome data. Can be: - Vector: named with sample IDs - Data
  frame: with row names matching bulk columns

- label_type:

  Character specifying phenotype label type (e.g., "SBS1", "time"),
  stored in `scRNA_data@misc`

- alpha:

  Parameter used to balance the effect of the l1 norm and the
  network-based penalties. It can be a number or a searching vector. If
  alpha = NULL, a default searching vector is used. The range of alpha
  is between 0 and 1. A larger alpha lays more emphasis on the l1 norm.

- cutoff:

  (default: `0.2`). When `alpha=NULL`, the cutoff is used to determine
  the optimal alpha. Higher values increase specificity.

- family:

  Model family for outcome type: - "gaussian": Continuous outcomes -
  "binomial": Binary outcomes (default) - "cox": Survival outcomes

- reliability_test:

  List controlling reliability testing:

  - run: Whether to perform reliability test (default: FALSE)

  - n: Permutation times (default: `10L`)

  - nfold: Cross-validation folds (default: `10L`)

- cell_evaluation:

  List controlling cell evaluation:

  - run: Whether to perform cell evaluation (default: FALSE)

  - benchmark_data: Path to benchmark data (RData file)

  - FDR_cutoff: FDR threshold for evaluation (default: `0.05`)

  - bootstrap_n: Bootstrap iterations (default: `100L`)

- ...:

  Additional arguments. Currently supports:

  - `verbose`: Logical indicating whether to print progress messages.
    Defaults to `TRUE`.

  - `seed`: For reproducibility, default is `123L`

  - `assay`: Assay to use for single-cell data. Defaults to \`"RNA"

## Value

A list containing:

- scRNA_data:

  A Seurat object with screened cells containing metadata:

  scissor

  :   "Positive"/"Negative"/"Neutral" classification

  label_type

  :   Outcome label used

- scissor_result:

  Raw Scissor results

- reliability_result:

  If reliability_test=TRUE, contains:

  statistic

  :   A value between 0 and 1

  p

  :   p-value of the test statistic

  AUC_test_real

  :   10 values of AUC for real data

  AUC_test_back

  :   A list of AUC for background data

- cell_evaluation:

  If cell_evaluation=TRUE, contains:

  evaluation_res

  :   A data.frame with some supporting information for each Scissor
      selected cell

## LICENSE

Licensed under the GNU General Public License version 3 (GPL-3.0). A
copy of the license is available at
<https://www.gnu.org/licenses/gpl-3.0.en.html>.

## References

Sun D, Guan X, Moran AE, Wu LY, Qian DZ, Schedin P, et al. Identifying
phenotype-associated subpopulations by integrating bulk and single-cell
sequencing data. Nat Biotechnol. 2022 Apr;40(4):527–38.

## See also

Other screen_method:
[`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md),
[`DoLP_SGL()`](https://wanglabcsu.github.io/sigbridger/reference/DoLP_SGL.md),
[`DoPIPET()`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md),
[`DoSIDISH()`](https://wanglabcsu.github.io/sigbridger/reference/DoSIDISH.md),
[`DoscAB()`](https://wanglabcsu.github.io/sigbridger/reference/DoscAB.md),
[`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md),
[`DoscPP()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPP.md)

Other scissor:
[`DoScissorCellEval()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissorCellEval.md),
[`DoScissorRelTest()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissorRelTest.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Binary outcome example
res <- DoScissor(
  matched_bulk = bulk_matrix,
  sc_data = seurat_obj,
  phenotype = a_named_vector,
  family = "binomial"
)
} # }
```
