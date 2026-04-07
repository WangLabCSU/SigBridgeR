# Perform PIPET Screening Analysis

Predicts cell subpopulations in single-cell data by matching expression
profiles to predefined marker gene templates using various
distance/similarity metrics. This function implements a template-based
classification approach with permutation testing for significance
assessment.

## Usage

``` r
DoPIPET(
  matched_bulk,
  sc_data,
  phenotype,
  phenotype_class = c("binary", "continuous", "survival"),
  only_pos_marker = FALSE,
  group = NULL,
  discretize_method = c("kmeans", "median", "custom"),
  cutoff = NULL,
  label_type = "PIPET",
  marker_finder = c("limma", "DESeq2"),
  log2FC = 1L,
  p_adjust = 0.05,
  show_log2FC = TRUE,
  freq_counts = NULL,
  normalize = TRUE,
  scale = TRUE,
  nPerm = 1000L,
  distance = c("cosine", "pearson", "spearman", "kendall", "euclidean", "maximum"),
  ...
)
```

## Arguments

- matched_bulk:

  Normalized bulk expression matrix (features × samples). Column names
  must match `phenotype` identifiers.

- sc_data:

  Seurat object containing single-cell RNA-seq data.

- phenotype:

  Clinical outcome data. Can be: - Vector: named with sample IDs - Data
  frame: with row names matching bulk columns

- phenotype_class:

  Analysis mode: - `"binary"`: Case-control design (e.g.,
  responder/non-responder) - `"continuous"`: Continuous outcome (e.g.,
  age, size) - `"survival"`: Patient survival

- only_pos_marker:

  A logical value. If `TRUE`, only upregulated marker genes (those with
  positive log2 fold change) will be retained. If `FALSE` (default),
  both upregulated and downregulated markers are kept.

- group:

  A character, name of one metadata column to group cells by (for
  example, orig.ident). The default value is `NULL`. In this case,
  screening will be performed on each group separately.

- discretize_method:

  `c("median", "kmeans", "custom")`. Discretization strategy for
  continuous phenotypes. Note: `"median"` is mapped internally to
  `"quantile"` (2-group quantile split). Default: `"kmeans"`.

- cutoff:

  Numeric vector of length `n_group - 1`. Required only when
  `discretize_method = "custom"`. Defines interior breakpoints on the
  *normalized, log2-transformed scale* (i.e., after
  `scale(log2(x + 1))`). Must be sorted in ascending order.

- label_type:

  Character specifying phenotype label type (e.g., "SBS1", "time"),
  stored in `scRNA_data@misc`

- marker_finder:

  A character, the marker finder method. The default value is `"limma"`.

- log2FC:

  In the DESeq differential expression analysis results, the cutoff
  value of log2FC. The default value is `1L`.

- p_adjust:

  In the DESeq differential expression analysis results, the cutoff
  value of adjust P. The default value is `0.05`.

- show_log2FC:

  Select whether to show log2 fold changes. The default value is `TRUE`.

- freq_counts:

  An integer, keep genes expressed in more than a certain number of
  cells. The default value is `NULL`, which means no filtering.

- normalize:

  Select whether to perform normalization of count data. The default
  value is `TRUE`.

- scale:

  Select whether to scale and center features in the dataset. The
  default value is `TRUE`.

- nPerm:

  An integer, number of permutations to do. The default value is
  `1000L`.

- distance:

  A character, the distance algorithm must be included in "cosine",
  "pearson", "spearman", "kendall","euclidean","maximum". default value
  is `NULL`, which means `"cosine"`.

- ...:

  Additional arguments to be passed to `PIPET.optimized`.

  - seed: Random seed for reproducibility

  - verbose: Whether to show progress messages

  - parallel: Whether to use parallel processing, default is `FALSE`.
    future::plan() must be set before calling this function.

  - assay: The assay to use, default is `"RNA"`

## See also

Other screen_method:
[`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md),
[`DoLP_SGL()`](https://wanglabcsu.github.io/sigbridger/reference/DoLP_SGL.md),
[`DoSCIPAC()`](https://wanglabcsu.github.io/sigbridger/reference/DoSCIPAC.md),
[`DoSIDISH()`](https://wanglabcsu.github.io/sigbridger/reference/DoSIDISH.md),
[`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md),
[`DoscAB()`](https://wanglabcsu.github.io/sigbridger/reference/DoscAB.md),
[`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md),
[`DoscPP()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPP.md)
