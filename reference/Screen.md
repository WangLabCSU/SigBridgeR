# Single-Cell Data Screening

Integrates matched bulk expression data and phenotype information to
identify phenotype-associated cell populations in single-cell RNA-seq
data using one of four computational methods. Ensures consistency
between bulk and phenotype data before analysis.

## Usage

``` r
Screen(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = NULL,
  phenotype_class = c("binary", "survival", "continuous"),
  screen_method = c("Scissor", "scPP", "scPAS", "scAB", "DEGAS", "LP_SGL", "PIPET"),
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
  columns), or - Patient survival Data frame with row names matching
  `matched_bulk` columns, colnames named "time" and "status"

- label_type:

  Character specifying phenotype label type (e.g., "SBS1", "time")

- phenotype_class:

  Type of phenotypic outcome (must be consistent with input data): -
  `"binary"`: Binary traits (e.g., case/control) - `"continuous"`:
  Continuous measurements - `"survival"`: Survival infomation

- screen_method:

  Screening algorithm to use, there are seven options:

  - "Scissor": see also DoScissor()

  - "scPP": see also DoscPP()

  - "scPAS": see also DoscPAS()

  - "scAB": see also DoscAB(), continuous phenotype is not supported

  - "DEGAS": see also DoDEGAS()

  - "LP_SGL": see also DoLP_SGL()

  - "PIPET": see also DoPIPET()

  - "SIDISH": see also DoSIDISH()

- ...:

  Additional method-specific parameters:

  Scissor

  :   

      alpha

      :   (numeric or NULL) Significance threshold. When NULL, alpha
          will keep increasing iteratively until the corresponding cells
          are screened out, default 0.05

      cutoff

      :   (numeric) A threshold for terminating the iteration of alpha,
          only work when alpha is NULL, default 0.2

      path2load_scissor_cache

      :   (character) default NULL

      path2save_scissor_inputs

      :   (character) A path to save the intermediary data. By using
          path2load_scissor_cache, the intermediary data can be loaded
          from the specified path. default "Scissor_inputs.RData"

      reliability_test

      :   (logical) Whether to perform reliability test, default FALSE

      reliability_test.nfold

      :   (integer) Cross-validation folds for reliability test, default
          10

      reliability_test.n

      :   (integer) Number of cells to use for reliability test, default
          10

      cell_evaluation

      :   (logical) Whether to perform cell evaluation, default FALSE

      cell_evaluation.benchmark_data

      :   .RData Benchmark data for cell evaluation, default NULL

      cell_evaluation.FDR

      :   (numeric) FDR threshold for cell evaluation, default 0.05

      cell_evaluation.bootstrap_n

      :   (integer) Number of bootstrap samples for cell evaluation,
          default 10

  scPP

  :   

      ref_group

      :   (integer or character) Reference group or baseline for binary
          comparisons, e.g. "Normal" for Tumor/Normal studies and 0 for
          0/1 case-control studies. default: 0

      Log2FC_cutoff

      :   (numeric) Minimum log2 fold-change for binary markers, default
          0.585

      estimate_cutoff

      :   (numeric) Effect size threshold for continuous traits, default
          0.2

      probs

      :   (numeric) Quantile cutoff for cell classification, default 0.2

  scPAS

  :   

      assay

      :   (character) Assay to use from sc_data, default "RNA"

      imputation

      :   (logical) Whether to perform imputation, default FALSE

      nfeature

      :   (integer) Number of features to select, default 3000

      alpha

      :   (numeric or NULL) Significance threshold, When NULL, alpha
          will keep increasing iteratively until the corresponding cells
          are screened out, default 0.01

      independent

      :   (logical) The background distribution of risk scores is
          constructed independently of each cell. default: TRUE

      network_class

      :   (character) Network class to use. default: 'SC', indicating
          gene-gene similarity networks derived from single-cell data.
          The other one is 'bulk'.

      permutation_times

      :   (integer) Number of permutations, default 2000

      FDR_threshold

      :   (numeric) FDR value threshold for identifying
          phenotype-associated cells default 0.05

  scAB

  :   

      alpha

      :   (numeric) Coefficient of phenotype regularization ,default
          0.005

      alpha_2

      :   (numeric) Coefficent of cell-cell similarity regularization,
          default 0.005

      maxiter

      :   (integer) NMF optimization iterations, default 2000

      tred

      :   (integer) Z-score threshold, default 2

  DEGAS

  :   

      sc_data.pheno_colname

      :   (character) Phenotype column name in sc_data, default "NULL"

      select_fraction

      :   (numeric) Fraction of cells to select for DEGAS, default 0.05

      tmp_dir

      :   (character) Temporary directory for DEGAS, default "NULL"

      env_params

      :   (list) Environment parameters for DEGAS, default "list()"

      degas_params

      :   (list) DEGAS parameters, default "list()"

      normality_test_method

      :   (character) Normality test method for DEGAS, default
          "jarque-bera"

  SIDISH

  :   

      sidish_params

      :   (list) SIDISH parameters, default "list()"

      env_params

      :   (list) Environment parameters for SIDISH, default "list()"

  LP_SGL

  :   

      resolution

      :   (numeric) Resolution parameter for Leiden clustering, default
          0.6

      alpha

      :   (numeric) Alpha parameter for SGL balancing L1 and L2
          penalties, default 0.5

      nfold

      :   (integer) Number of folds for cross-validation, default 5

      dge_analysis

      :   (list) Differential expression analysis settings:

          - run: (logical) Whether to run DEG analysis, default FALSE

          - logFC_threshold: (numeric) Log fold change threshold,
            default 1

          - pval_threshold: (numeric) P-value threshold, default 0.05

      PIPET

      :   

          group

          :   (character or NULL) Name of a metadata column (e.g.,
              `"orig.ident"`) to stratify cells before screening. When
              `NULL` (default), screening is performed globally across
              all cells.

          discretize_method

          :   (character) Strategy to binarize continuous phenotypes
              internally before marker identification. One of:

              - `"median"` (default): Equivalent to 2-quantile split
                (i.e., median threshold).

              - `"kmeans"`: Two-cluster k-means on the continuous
                phenotype.

              - `"custom"`: User-defined cutoffs via `cutoff`.

          cutoff

          :   (numeric vector or NULL) Required only if
              `discretize_method = "custom"`. Specifies interior
              breakpoints on the *normalized, log2-transformed phenotype
              scale* (i.e., after `scale(log2(x + 1))`). Must be sorted
              ascending and of length `n_group - 1`.

          label_type

          :   (character) Phenotype label type (e.g., `"PIPET_SBS1"`),
              stored in `scRNA_data@misc`. Default: `"PIPET"`.

          log2FC

          :   (numeric) Absolute log2 fold-change cutoff for
              differential expression marker selection in bulk data (via
              DESeq2-like analysis). Default: 1.

          p_adjust

          :   (numeric) Adjusted p-value (FDR) cutoff for marker gene
              selection. Default: 0.05.

          show_log2FC

          :   (logical) Whether to annotate markers with signed log2FC
              direction (e.g., `CD3D_up`). Default: `TRUE`.

          freq_counts

          :   (integer or NULL) Minimum number of cells a gene must be
              expressed in to be retained in scRNA-seq data
              preprocessing. Default: `NULL` (no filtering).

          normalize

          :   (logical) Whether to apply log-normalization
              (`LogNormalize`) to scRNA-seq counts prior to correlation.
              Default: `TRUE`.

          scale

          :   (logical) Whether to scale (center + unit-variance) gene
              expression across cells before computing distances.
              Default: `TRUE`.

          nPerm

          :   (integer) Number of label permutations to assess
              significance of correlation scores. Default: 1000.

          distance

          :   (character) Distance or similarity metric for template
              matching. Supported: `"cosine"` (default), `"pearson"`,
              `"spearman"`, `"kendall"`, `"euclidean"`, `"maximum"`.

          seed

          :   (integer or NULL) Random seed for reproducibility in
              marker creation and permutation tests. Default: inherits
              from `getFuncOption("seed")`.

          verbose

          :   (logical) Whether to print progress messages. Default:
              inherits from `getFuncOption("verbose")`.

          parallel

          :   (logical) Whether to enable parallel permutations
              (requires
              [`future::plan()`](https://future.futureverse.org/reference/plan.html)
              pre-set). Default: `FALSE`.

## Value

A list containing:

- scRNA_data:

  A Seurat object with phenotype-associated cells labelled in
  `meta.data` column

- Some screen_result:

  Important information about the screened result related to the
  selected method

## Data Matching Requirements

- matched_bulk column names and phenotype names/rownames must be
  identical

- Phenotype values must correspond to bulk samples (not directly to
  single cells)

- Mismatches will trigger an error before analysis begins, and there is
  a built-in pre-run check.

## Method Compatibility

|         |                      |                                                                                                                                                                                                                               |
|---------|----------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Method  | Supported Phenotypes | Additional Parameters                                                                                                                                                                                                         |
| Scissor | All three types      | alpha, cutoff, path2load_scissor_cache, path2save_scissor_inputs, reliability_test, reliability_test.n,reliability_test.nfold, cell_evaluation,cell_evaluation.benchmark_data,cell_evaluation.FDR,cell_evaluation.bootstrap_n |
| scPP    | All three types      | ref_group, Log2FC_cutoff, estimate_cutoff, probs                                                                                                                                                                              |
| scPAS   | All three types      | n_components ,assay, imputation,nfeature, alpha,network_class,permutation_times,FDR_threshold,independent                                                                                                                     |
| scAB    | Binary/Survival      | alpha, alpha_2, maxiter, tred                                                                                                                                                                                                 |
| DEGAS   | All three types      | sc_data.pheno_colname,select_fraction,tmp_dir,env_params,degas_params,normality_test_method                                                                                                                                   |
| LP_SGL  | All three types      | resolution, alpha, nfold, dge_analysis                                                                                                                                                                                        |
| PIPET   | Binary/Continuous    | group, discretize_method, cutoff, log2FC, p_adjust, show_log2FC, freq_counts, normalize, scale, nPerm, distance                                                                                                               |
| SIDISH  | Survival Only        | sidish_params, env_params                                                                                                                                                                                                     |

## See also

Associated functions:

- [`DoScissor`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md)

- [`DoscPP`](https://wanglabcsu.github.io/sigbridger/reference/DoscPP.md)

- [`DoscPAS`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md)

- [`DoscAB`](https://wanglabcsu.github.io/sigbridger/reference/DoscAB.md)

- [`DoDEGAS`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md)

- [`DoLP_SGL`](https://wanglabcsu.github.io/sigbridger/reference/DoLP_SGL.md)

- [`DoPIPET`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md)

- [`DoSIDISH`](https://wanglabcsu.github.io/sigbridger/reference/DoSIDISH.md)
