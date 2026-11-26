# Single-Cell Data Screening

Integrates matched bulk expression data and phenotype information to
identify phenotype-associated cell populations in single-cell RNA-seq
data using a registered screening method. Ensures consistency between
bulk and phenotype data before analysis.

## Usage

``` r
Screen(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = NULL,
  phenotype_class = c("binary", "survival", "continuous"),
  screen_method = c("Scissor", "scPP", "scPAS", "scAB", "DEGAS", "LP_SGL", "PIPET",
    "SIDISH", "SCIPAC", "TiRank", "DEGASv2", "PhenoMapR"),
  ...
)
```

## Arguments

- matched_bulk:

  Matrix or data frame of preprocessed bulk RNA-seq expression data
  (genes x samples). Column names must match names/IDs in phenotype.

- sc_data:

  A Seurat object containing scRNA-seq data to be screened.

- phenotype:

  Phenotype data, either:

  - Named vector, where names match columns of matched_bulk

  - Patient survival data frame, where row names match columns of
    matched_bulk, and column names are "time" and "status"

- label_type:

  Character specifying phenotype label type, e.g. "SBS1" or "time".

- phenotype_class:

  Type of phenotypic outcome. Must be consistent with input data:

  - "binary": Binary traits, e.g. case/control

  - "continuous": Continuous measurements

  - "survival": Survival information

- screen_method:

  Screening algorithm to use. Available options:

  - "DEGASv2": see also
    [`DoDEGASv2()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGASv2.md)

  - "Scissor": see also
    [`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md)

  - "TiRank": see also
    [`DoTiRank()`](https://wanglabcsu.github.io/sigbridger/reference/DoTiRank.md)

  - "SCIPAC": see also
    [`DoSCIPAC()`](https://wanglabcsu.github.io/sigbridger/reference/DoSCIPAC.md)

  - "scPP": see also
    [`DoscPP()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPP.md)

  - "LP_SGL": see also
    [`DoLP_SGL()`](https://wanglabcsu.github.io/sigbridger/reference/DoLP_SGL.md)

  - "PIPET": see also
    [`DoPIPET()`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md)

  - "DEGAS": see also
    [`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md)

  - "scAB": see also
    [`DoscAB()`](https://wanglabcsu.github.io/sigbridger/reference/DoscAB.md)

  - "PhenoMapR": see also
    [`DoPhenoMapR()`](https://wanglabcsu.github.io/sigbridger/reference/DoPhenoMapR.md)

  - "SIDISH": see also
    [`DoSIDISH()`](https://wanglabcsu.github.io/sigbridger/reference/DoSIDISH.md)

  - "scPAS": see also
    [`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md)

- ...:

  Additional method-specific parameters.

  Scissor

  - alpha: Numeric or NULL. Significance threshold. When NULL, alpha
    increases iteratively until corresponding cells are screened out.
    Default: 0.05.

  - cutoff: Numeric. Threshold for terminating the iteration of alpha;
    only works when alpha = NULL. Default: 0.2.

  - path2load_scissor_cache: Character. Default: NULL.

  - path2save_scissor_inputs: Character. Path to save intermediary data.
    Can be loaded later via path2load_scissor_cache. Default:
    "Scissor_inputs.RData".

  - reliability_test: Logical. Whether to perform reliability test.
    Default: FALSE.

  - reliability_test.nfold: Integer. Cross-validation folds for
    reliability test. Default: 10.

  - reliability_test.n: Integer. Number of cells to use for reliability
    test. Default: 10.

  - cell_evaluation: Logical. Whether to perform cell evaluation.
    Default: FALSE.

  - cell_evaluation.benchmark_data: .RData benchmark data for cell
    evaluation. Default: NULL.

  - cell_evaluation.FDR: Numeric. FDR threshold for cell evaluation.
    Default: 0.05.

  - cell_evaluation.bootstrap_n: Integer. Number of bootstrap samples
    for cell evaluation. Default: 10.

  scPP

  - ref_group: Integer or character. Reference group or baseline for
    binary comparisons, e.g. "Normal" or 0. Default: 0.

  - Log2FC_cutoff: Numeric. Minimum log2 fold-change for binary markers.
    Default: 0.585.

  - estimate_cutoff: Numeric. Effect size threshold for continuous
    traits. Default: 0.2.

  - probs: Numeric. Quantile cutoff for cell classification. Default:
    0.2.

  scPAS

  - assay: Character. Assay to use from sc_data. Default: "RNA".

  - imputation: Logical. Whether to perform imputation. Default: FALSE.

  - nfeature: Integer. Number of features to select. Default: 3000.

  - alpha: Numeric or NULL. Significance threshold. When NULL, alpha
    increases iteratively until corresponding cells are screened out.
    Default: 0.01.

  - independent: Logical. Whether the background distribution of risk
    scores is constructed independently of each cell. Default: TRUE.

  - network_class: Character. Network class to use. Default: "SC" for
    single-cell-derived gene-gene similarity networks; another option is
    "bulk".

  - permutation_times: Integer. Number of permutations. Default: 2000.

  - FDR_threshold: Numeric. FDR threshold for identifying
    phenotype-associated cells. Default: 0.05.

  scAB

  - alpha: Numeric. Coefficient of phenotype regularization. Default:
    0.005.

  - alpha_2: Numeric. Coefficient of cell-cell similarity
    regularization. Default: 0.005.

  - maxiter: Integer. NMF optimization iterations. Default: 2000.

  - tred: Integer. Z-score threshold. Default: 2.

  DEGAS

  - sc_data.pheno_colname: Character. Phenotype column name in sc_data.
    Default: "NULL".

  - select_fraction: Numeric. Fraction of cells to select for DEGAS.
    Default: 0.05.

  - tmp_dir: Character. Temporary directory for DEGAS. Default: "NULL".

  - env_params: List. Environment parameters for DEGAS. Default: list().

  - degas_params: List. DEGAS parameters. Default: list().

  - normality_test_method: Character. Normality test method for DEGAS.
    Default: "jarque-bera".

  SIDISH

  - sidish_params: List. SIDISH parameters. Default: list().

  - env_params: List. Environment parameters for SIDISH. Default:
    list().

  SCIPAC

  - hvg: Integer. Number of highly variable genes to use for
    preprocessing. Default: 1000.

  - do_pca_sc: Logical. If TRUE, first do PCA on sc.dat and use the
    rotation matrix on bulk.dat; if FALSE, first do PCA on bulk.data and
    use the rotation matrix on sc.dat. Default: FALSE.

  - n_pc: Integer. Number of principal components to use. Default: 60.

  - sc_batch_col: Character or vector. Batch variable for single-cell
    data. If character, should be a column name in `sc_data@metadata`.
    If vector, should match cell order. Default: NULL.

  - resolution: Integer. Clustering resolution for cell type
    identification. Higher values produce more clusters. Default: 2.

  - ela_net_alpha: Numeric. Elastic net mixing parameter (0 = ridge, 1 =
    lasso). Default: 0.4.

  - bt_size: Integer. Bootstrap sample size for stability assessment.
    Default: 50.

  - ncore: Integer. Number of CPU cores for parallel computation.
    Default: 7.

  - ci_alpha: Numeric. Significance level for confidence intervals.
    Default: 0.05.

  - nfold: Integer. Number of folds for cross-validation. Default: 10.

  - assay: Character. Assay to use from sc_data. Default: "RNA".

  - verbose: Logical. Whether to print progress messages. Default:
    inherits from getFuncOption("verbose").

  - seed: Integer. Random seed for reproducibility. Default: 123.

  PIPET

  - group: Character or NULL. Metadata column used to stratify cells
    before screening. Default: NULL.

  - discretize_method: Character. Strategy to binarize continuous
    phenotypes:

    - "median": Median threshold. Default.

    - "kmeans": Two-cluster k-means.

    - "custom": User-defined cutoffs via cutoff.

  - cutoff: Numeric vector or NULL. Required when discretize_method =
    "custom". Specifies interior breakpoints on the normalized,
    log2-transformed phenotype scale, after scale(log2(x + 1)).

  - label_type: Character. Phenotype label type stored in
    scRNA_data@misc. Default: "PIPET".

  - log2FC: Numeric. Absolute log2 fold-change cutoff for bulk marker
    selection. Default: 1.

  - p_adjust: Numeric. Adjusted p-value/FDR cutoff for marker selection.
    Default: 0.05.

  - show_log2FC: Logical. Whether to annotate markers with signed log2FC
    direction. Default: TRUE.

  - freq_counts: Integer or NULL. Minimum number of cells a gene must be
    expressed in. Default: NULL.

  - normalize: Logical. Whether to apply LogNormalize before
    correlation. Default: TRUE.

  - scale: Logical. Whether to center and scale gene expression before
    computing distances. Default: TRUE.

  - nPerm: Integer. Number of label permutations. Default: 1000.

  - distance: Character. Supported metrics: "cosine", "pearson",
    "spearman", "kendall", "euclidean", and "maximum". Default:
    "cosine".

  - seed: Integer or NULL. Random seed. Default: inherits from
    getFuncOption("seed").

  - verbose: Logical. Whether to print progress messages. Default:
    inherits from getFuncOption("verbose").

  - parallel: Logical. Whether to enable parallel permutations. Requires
    future::plan() to be pre-set. Default: FALSE.

  LP_SGL

  - resolution: Numeric. Resolution parameter for Leiden clustering.
    Default: 0.6.

  - alpha: Numeric. SGL parameter balancing L1 and L2 penalties.
    Default: 0.5.

  - nfold: Integer. Number of folds for cross-validation. Default: 5.

  - dge_analysis: List. Differential expression analysis settings:

    - run: Logical. Whether to run DEG analysis. Default: FALSE.

    - logFC_threshold: Numeric. Log fold-change threshold. Default: 1.

    - pval_threshold: Numeric. P-value threshold. Default: 0.05.

  DEGASv2

  - top_fraction_pos: Numeric. Proportion of cells to be labelled as
    "Positive". Default: 0.2.

  - sclab: Character or NULL. Metadata column providing supervised cell
    labels. Default: NULL.

  - bulk_hvg: Logical. Whether to select highly variable genes from bulk
    data. Default: TRUE.

  - bulk_de: Logical. Whether to perform differential expression
    analysis on bulk data. Default: TRUE.

  - sc_de: Logical. Whether to perform differential expression analysis
    on single-cell data. Default: TRUE.

  - add_genes: Character vector or NULL. Additional genes to include in
    the analysis. Default: NULL.

  - n_hvg: Integer. Number of highly variable genes to select. Default:
    250.

  - n_bulk_de: Integer. Number of bulk DE genes to select. Default: 250.

  - n_sc_de: Integer. Number of single-cell DE genes to select. Default:
    200.

  - padj.thresh: Numeric. Adjusted p-value threshold for DE gene
    selection. Default: 0.05.

  - only.pos: Logical. Whether to keep only positively enriched markers.
    Default: FALSE.

  - min.pct: Numeric. Minimum fraction of cells expressing a gene for DE
    analysis. Default: 0.25.

  - logfc.threshold: Numeric. Log fold-change threshold for DE analysis.
    Default: 0.25.

  - n_st_classes: Integer. Number of single-cell classes. Default:
    length(unique(sc_data\$seurat_clusters)).

  - loss_type: Character. Loss function for model training. One of
    "cross_entropy", "log_neg", "rank_loss". Default: "cross_entropy".

  - transfer_type: Character. Domain adaptation method. One of
    "Wasserstein", "MMD". Default: "Wasserstein".

  - model_save_dir: Character. Directory to save the trained model.
    Default: "DEGASv2_res".

  - lambda1: Numeric. Weight of the classification loss. Default: 1.

  - lambda2: Numeric. Weight of the domain-adaptation loss. Default: 3.

  - lambda3: Numeric. Weight of the reconstruction loss. Default: 3.

  - tot_seeds: Integer. Number of random seeds to train. Default: 10.

  - tot_iters: Integer. Number of training iterations. Default: 300.

  - extract_embs: Logical. Whether to extract cell embeddings. Default:
    FALSE.

  - random_feat: Logical. Whether to use random features. Default:
    FALSE.

  - random_perc: Numeric. Proportion of random features when random_feat
    = TRUE. Default: 0.8.

  - early_stopping: Logical. Whether to enable early stopping. Default:
    FALSE.

  TiRank

  - tirank_params: List. TiRank algorithm parameters, including data
    preprocessing (validation_proportion, sampling_thresh,
    sampling_mode, top_var_genes, top_gene_pairs, p_value_threshold,
    max_cutoff, min_cutoff), neural network architecture (nhead, nhid1,
    nhid2, n_output, nlayers, n_pred, dropout, encoder_type,
    infer_mode), and training (n_trials, do_reject, tolerance). See
    [`DoTiRank()`](https://wanglabcsu.github.io/sigbridger/reference/DoTiRank.md)
    for details. Default: list().

  - save_path: Character. Directory to save intermediate and final
    results. Soft-deprecated; prefer passing load_cache/save_cache via
    ... Default: "./TiRank_res".

  - load_cache: Character or NULL. Path to cached TiRank data.
    Soft-deprecated; prefer passing load_cache via ... Default: NULL.

  - verbose: Logical. Whether to print progress messages. Default:
    inherits from getFuncOption("verbose").

  - seed: Integer. Random seed for reproducibility. Default: inherits
    from getFuncOption("seed").

  - assay: Character. Assay to use from sc_data. Default: "RNA".

  - save_cache: Character. Directory to save TiRank cache. Default:
    save_path.

  - additional_description: Character. Optional description appended to
    the cache metadata.

  PhenoMapR

  - phenotype_column: Character or NULL. Outcome column name for
    binary/continuous phenotypes. Default: NULL (uses
    "phenotype_value").

  - survival_time, survival_event: Character. Columns of the survival
    phenotype data frame. Default: "time" / "status".

  - normalize: Logical. Whether to normalize bulk expression. Default:
    TRUE.

  - platform: Character. "auto", "rnaseq", or "microarray". Default:
    "auto".

  - probe_annotation: Probe-to-gene annotation for microarray data.
    Default: NULL.

  - hugo_species: Character vector passed to HGNChelper. Default:
    c("human", "mouse").

  - binary_positive_reference: Character. "second" or "first". Default:
    "second".

  - assay, layer: Character. Seurat assay and layer used for scoring.
    Default: "RNA" / "data".

  - z_score_cutoff: Numeric. Z-score cutoff for limiting paired genes.
    Default: 2.

  - reference_sign: Integer. Sign applied to reference scores. Default:
    1.

  - score_mode: Character vector. Default: c("weighted_sum",
    "activity_adjusted").

  - vars_to_regress: Character vector. Default: c("S.Score",
    "G2M.Score", "nCount_RNA").

  - cancer_type: Character or NULL. Optional cancer-type context.
    Default: NULL.

  - pseudobulk: Logical. Whether to score pseudobulk profiles. Default:
    FALSE.

  - group_by: Character or NULL. Grouping variable for pseudobulk.
    Default: NULL.

  - sample_metadata: Data frame or NULL. Metadata for pseudobulk.
    Default: NULL.

  - percentile: Numeric. Tail percentile for phenotype groups. Default:
    0.05.

  - permutation_n: Integer or NULL. Number of permutations. Default:
    NULL.

  - permutation_seed: Integer or NULL. Permutation seed; defaults to
    seed.

## Value

A list containing:

- scRNA_data: A Seurat object with phenotype-associated cells labelled
  in a meta.data column.

- screen_result: Important information about the screened result related
  to the selected method.

## Data Matching Requirements

- Column names of matched_bulk and names or row names of phenotype must
  be identical.

- Phenotype values must correspond to bulk samples, not directly to
  single cells.

- Mismatches will trigger an error before analysis begins.

- A built-in pre-run check is performed.

## See also

Associated functions:

- [`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md)

- [`DoscPP()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPP.md)

- [`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md)

- [`DoscAB()`](https://wanglabcsu.github.io/sigbridger/reference/DoscAB.md)

- [`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md)

- [`DoLP_SGL()`](https://wanglabcsu.github.io/sigbridger/reference/DoLP_SGL.md)

- [`DoPIPET()`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md)

- [`DoSIDISH()`](https://wanglabcsu.github.io/sigbridger/reference/DoSIDISH.md)

- [`DoSCIPAC()`](https://wanglabcsu.github.io/sigbridger/reference/DoSCIPAC.md)

- [`DoTiRank()`](https://wanglabcsu.github.io/sigbridger/reference/DoTiRank.md)

- [`DoDEGASv2()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGASv2.md)

- [`DoPhenoMapR()`](https://wanglabcsu.github.io/sigbridger/reference/DoPhenoMapR.md)
