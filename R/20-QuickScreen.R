#' @title Single-Cell Data Screening
#'
#' @description
#' Integrates matched bulk expression data and phenotype information to identify
#' phenotype-associated cell populations in single-cell RNA-seq data using one of
#' nine computational methods. Ensures consistency between bulk and phenotype data
#' before analysis.
#'
#' @param matched_bulk Matrix or data frame of preprocessed bulk RNA-seq expression
#' data (genes x samples). Column names must match names/IDs in phenotype.
#' @param sc_data A Seurat object containing scRNA-seq data to be screened.
#' @param phenotype Phenotype data, either:
#'
#' - Named vector, where names match columns of matched_bulk
#' - Patient survival data frame, where row names match columns of matched_bulk,
#'   and column names are "time" and "status"
#'
#' @param label_type Character specifying phenotype label type, e.g. "SBS1" or
#' "time".
#' @param phenotype_class Type of phenotypic outcome. Must be consistent with input
#' data:
#'
#' - "binary": Binary traits, e.g. case/control
#' - "continuous": Continuous measurements
#' - "survival": Survival information
#'
#' @param screen_method Screening algorithm to use. Available options:
#'
#' - "Scissor": see also [DoScissor()]
#' - "scPP": see also [DoscPP()]
#' - "scPAS": see also [DoscPAS()]
#' - "scAB": see also [DoscAB()]; continuous phenotype is not supported
#' - "DEGAS": see also [DoDEGAS()]
#' - "DEGASv2": see also [DoDEGASv2()]
#' - "LP_SGL": see also [DoLP_SGL()]
#' - "PIPET": see also [DoPIPET()]
#' - "SIDISH": see also [DoSIDISH()]
#' - "SCIPAC": see also [DoSCIPAC()]
#' - "TiRank": see also [DoTiRank()]
#'
#' @param ... Additional method-specific parameters.
#'
#' Scissor
#'
#' - alpha: Numeric or NULL. Significance threshold. When NULL, alpha
#'   increases iteratively until corresponding cells are screened out. Default: 0.05.
#' - cutoff: Numeric. Threshold for terminating the iteration of alpha; only
#'   works when alpha = NULL. Default: 0.2.
#' - path2load_scissor_cache: Character. Default: NULL.
#' - path2save_scissor_inputs: Character. Path to save intermediary data. Can be
#'   loaded later via path2load_scissor_cache. Default: "Scissor_inputs.RData".
#' - reliability_test: Logical. Whether to perform reliability test. Default: FALSE.
#' - reliability_test.nfold: Integer. Cross-validation folds for reliability test.
#'   Default: 10.
#' - reliability_test.n: Integer. Number of cells to use for reliability test.
#'   Default: 10.
#' - cell_evaluation: Logical. Whether to perform cell evaluation. Default: FALSE.
#' - cell_evaluation.benchmark_data: .RData benchmark data for cell evaluation.
#'   Default: NULL.
#' - cell_evaluation.FDR: Numeric. FDR threshold for cell evaluation. Default: 0.05.
#' - cell_evaluation.bootstrap_n: Integer. Number of bootstrap samples for cell
#'   evaluation. Default: 10.
#'
#' scPP
#'
#' - ref_group: Integer or character. Reference group or baseline for binary
#'   comparisons, e.g. "Normal" or 0. Default: 0.
#' - Log2FC_cutoff: Numeric. Minimum log2 fold-change for binary markers.
#'   Default: 0.585.
#' - estimate_cutoff: Numeric. Effect size threshold for continuous traits.
#'   Default: 0.2.
#' - probs: Numeric. Quantile cutoff for cell classification. Default: 0.2.
#'
#' scPAS
#'
#' - assay: Character. Assay to use from sc_data. Default: "RNA".
#' - imputation: Logical. Whether to perform imputation. Default: FALSE.
#' - nfeature: Integer. Number of features to select. Default: 3000.
#' - alpha: Numeric or NULL. Significance threshold. When NULL, alpha
#'   increases iteratively until corresponding cells are screened out. Default: 0.01.
#' - independent: Logical. Whether the background distribution of risk scores is
#'   constructed independently of each cell. Default: TRUE.
#' - network_class: Character. Network class to use. Default: "SC" for
#'   single-cell-derived gene-gene similarity networks; another option is "bulk".
#' - permutation_times: Integer. Number of permutations. Default: 2000.
#' - FDR_threshold: Numeric. FDR threshold for identifying phenotype-associated
#'   cells. Default: 0.05.
#'
#' scAB
#'
#' - alpha: Numeric. Coefficient of phenotype regularization. Default: 0.005.
#' - alpha_2: Numeric. Coefficient of cell-cell similarity regularization.
#'   Default: 0.005.
#' - maxiter: Integer. NMF optimization iterations. Default: 2000.
#' - tred: Integer. Z-score threshold. Default: 2.
#'
#' DEGAS
#'
#' - sc_data.pheno_colname: Character. Phenotype column name in sc_data.
#'   Default: "NULL".
#' - select_fraction: Numeric. Fraction of cells to select for DEGAS. Default: 0.05.
#' - tmp_dir: Character. Temporary directory for DEGAS. Default: "NULL".
#' - env_params: List. Environment parameters for DEGAS. Default: list().
#' - degas_params: List. DEGAS parameters. Default: list().
#' - normality_test_method: Character. Normality test method for DEGAS.
#'   Default: "jarque-bera".
#'
#' SIDISH
#'
#' - sidish_params: List. SIDISH parameters. Default: list().
#' - env_params: List. Environment parameters for SIDISH. Default: list().
#'
#' SCIPAC
#'
#' - hvg: Integer. Number of highly variable genes to use for preprocessing.
#'   Default: 1000.
#' - do_pca_sc: Logical. If TRUE, first do PCA on sc.dat and use the rotation
#'   matrix on bulk.dat; if FALSE, first do PCA on bulk.data and use the
#'   rotation matrix on sc.dat. Default: FALSE.
#' - n_pc: Integer. Number of principal components to use. Default: 60.
#' - sc_batch_col: Character or vector. Batch variable for single-cell data.
#'   If character, should be a column name in `sc_data@metadata`. If vector, should
#'   match cell order. Default: NULL.
#' - resolution: Integer. Clustering resolution for cell type identification.
#'   Higher values produce more clusters. Default: 2.
#' - ela_net_alpha: Numeric. Elastic net mixing parameter (0 = ridge,
#'   1 = lasso). Default: 0.4.
#' - bt_size: Integer. Bootstrap sample size for stability assessment. Default: 50.
#' - ncore: Integer. Number of CPU cores for parallel computation. Default: 7.
#' - ci_alpha: Numeric. Significance level for confidence intervals. Default: 0.05.
#' - nfold: Integer. Number of folds for cross-validation. Default: 10.
#' - assay: Character. Assay to use from sc_data. Default: "RNA".
#' - verbose: Logical. Whether to print progress messages. Default: inherits from
#'   getFuncOption("verbose").
#' - seed: Integer. Random seed for reproducibility. Default: 123.
#'
#' PIPET
#'
#' - group: Character or NULL. Metadata column used to stratify cells before
#'   screening. Default: NULL.
#' - discretize_method: Character. Strategy to binarize continuous phenotypes:
#'   - "median": Median threshold. Default.
#'   - "kmeans": Two-cluster k-means.
#'   - "custom": User-defined cutoffs via cutoff.
#' - cutoff: Numeric vector or NULL. Required when discretize_method = "custom".
#'   Specifies interior breakpoints on the normalized, log2-transformed phenotype
#'   scale, after scale(log2(x + 1)).
#' - label_type: Character. Phenotype label type stored in scRNA_data@misc.
#'   Default: "PIPET".
#' - log2FC: Numeric. Absolute log2 fold-change cutoff for bulk marker selection.
#'   Default: 1.
#' - p_adjust: Numeric. Adjusted p-value/FDR cutoff for marker selection.
#'   Default: 0.05.
#' - show_log2FC: Logical. Whether to annotate markers with signed log2FC direction.
#'   Default: TRUE.
#' - freq_counts: Integer or NULL. Minimum number of cells a gene must be
#'   expressed in. Default: NULL.
#' - normalize: Logical. Whether to apply LogNormalize before correlation.
#'   Default: TRUE.
#' - scale: Logical. Whether to center and scale gene expression before computing
#'   distances. Default: TRUE.
#' - nPerm: Integer. Number of label permutations. Default: 1000.
#' - distance: Character. Supported metrics: "cosine", "pearson", "spearman",
#'   "kendall", "euclidean", and "maximum". Default: "cosine".
#' - seed: Integer or NULL. Random seed. Default: inherits from
#'   getFuncOption("seed").
#' - verbose: Logical. Whether to print progress messages. Default: inherits from
#'   getFuncOption("verbose").
#' - parallel: Logical. Whether to enable parallel permutations. Requires
#'   future::plan() to be pre-set. Default: FALSE.
#'
#' LP_SGL
#'
#' - resolution: Numeric. Resolution parameter for Leiden clustering. Default: 0.6.
#' - alpha: Numeric. SGL parameter balancing L1 and L2 penalties. Default: 0.5.
#' - nfold: Integer. Number of folds for cross-validation. Default: 5.
#' - dge_analysis: List. Differential expression analysis settings:
#'   - run: Logical. Whether to run DEG analysis. Default: FALSE.
#'   - logFC_threshold: Numeric. Log fold-change threshold. Default: 1.
#'   - pval_threshold: Numeric. P-value threshold. Default: 0.05.
#'
#' DEGASv2
#'
#' - top_fraction_pos: Numeric. Proportion of cells to be labelled as
#'   "Positive". Default: 0.2.
#' - sclab: Character or NULL. Metadata column providing supervised cell
#'   labels. Default: NULL.
#' - bulk_hvg: Logical. Whether to select highly variable genes from bulk
#'   data. Default: TRUE.
#' - bulk_de: Logical. Whether to perform differential expression analysis on
#'   bulk data. Default: TRUE.
#' - sc_de: Logical. Whether to perform differential expression analysis on
#'   single-cell data. Default: TRUE.
#' - add_genes: Character vector or NULL. Additional genes to include in the
#'   analysis. Default: NULL.
#' - n_hvg: Integer. Number of highly variable genes to select. Default: 250.
#' - n_bulk_de: Integer. Number of bulk DE genes to select. Default: 250.
#' - n_sc_de: Integer. Number of single-cell DE genes to select. Default: 200.
#' - padj.thresh: Numeric. Adjusted p-value threshold for DE gene selection.
#'   Default: 0.05.
#' - only.pos: Logical. Whether to keep only positively enriched markers.
#'   Default: FALSE.
#' - min.pct: Numeric. Minimum fraction of cells expressing a gene for DE
#'   analysis. Default: 0.25.
#' - logfc.threshold: Numeric. Log fold-change threshold for DE analysis.
#'   Default: 0.25.
#' - n_st_classes: Integer. Number of single-cell classes. Default:
#'   length(unique(sc_data$seurat_clusters)).
#' - loss_type: Character. Loss function for model training. One of
#'   "cross_entropy", "log_neg", "rank_loss". Default: "cross_entropy".
#' - transfer_type: Character. Domain adaptation method. One of "Wasserstein",
#'   "MMD". Default: "Wasserstein".
#' - model_save_dir: Character. Directory to save the trained model. Default:
#'   "DEGASv2_res".
#' - lambda1: Numeric. Weight of the classification loss. Default: 1.
#' - lambda2: Numeric. Weight of the domain-adaptation loss. Default: 3.
#' - lambda3: Numeric. Weight of the reconstruction loss. Default: 3.
#' - tot_seeds: Integer. Number of random seeds to train. Default: 10.
#' - tot_iters: Integer. Number of training iterations. Default: 300.
#' - extract_embs: Logical. Whether to extract cell embeddings. Default: FALSE.
#' - random_feat: Logical. Whether to use random features. Default: FALSE.
#' - random_perc: Numeric. Proportion of random features when random_feat =
#'   TRUE. Default: 0.8.
#' - early_stopping: Logical. Whether to enable early stopping. Default: FALSE.
#'
#' TiRank
#'
#' - tirank_params: List. TiRank algorithm parameters, including data
#'   preprocessing (validation_proportion, sampling_thresh, sampling_mode,
#'   top_var_genes, top_gene_pairs, p_value_threshold, max_cutoff, min_cutoff),
#'   neural network architecture (nhead, nhid1, nhid2, n_output, nlayers,
#'   n_pred, dropout, encoder_type, infer_mode), and training (n_trials,
#'   do_reject, tolerance). See [DoTiRank()] for details. Default: list().
#' - save_path: Character. Directory to save intermediate and final results.
#'   Soft-deprecated; prefer passing load_cache/save_cache via ...
#'   Default: "./TiRank_res".
#' - load_cache: Character or NULL. Path to cached TiRank data. Soft-deprecated;
#'   prefer passing load_cache via ... Default: NULL.
#' - verbose: Logical. Whether to print progress messages. Default: inherits
#'   from getFuncOption("verbose").
#' - seed: Integer. Random seed for reproducibility. Default: inherits from
#'   getFuncOption("seed").
#' - assay: Character. Assay to use from sc_data. Default: "RNA".
#' - save_cache: Character. Directory to save TiRank cache. Default: save_path.
#' - additional_description: Character. Optional description appended to the
#'   cache metadata.
#'
#' @return A list containing:
#'
#' - scRNA_data: A Seurat object with phenotype-associated cells labelled in a
#'   meta.data column.
#' - screen_result: Important information about the screened result related to
#'   the selected method.
#'
#' @section Data Matching Requirements:
#'
#' - Column names of matched_bulk and names or row names of phenotype must be
#'   identical.
#' - Phenotype values must correspond to bulk samples, not directly to single cells.
#' - Mismatches will trigger an error before analysis begins.
#' - A built-in pre-run check is performed.
#'
#'
#' @seealso
#' Associated functions:
#'
#' - [DoScissor()]
#' - [DoscPP()]
#' - [DoscPAS()]
#' - [DoscAB()]
#' - [DoDEGAS()]
#' - [DoLP_SGL()]
#' - [DoPIPET()]
#' - [DoSIDISH()]
#' - [DoSCIPAC()]
#' - [DoTiRank()]
#' - [DoDEGASv2()]
#'
#' @export
Screen <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = NULL,
  phenotype_class = c("binary", "survival", "continuous"),
  screen_method = c(
    "Scissor",
    "scPP",
    "scPAS",
    "scAB",
    "DEGAS",
    "LP_SGL",
    "PIPET",
    "SIDISH",
    "SCIPAC",
    "TiRank",
    "DEGASv2"
  ),
  ...
) {
  on.exit(gc(verbose = FALSE))

  if (is.null(label_type) || length(label_type) != 1) {
    cli::cli_alert_info(c(
      "i" = "{.var label_type} not specified or not of length 1, using {.val {screen_method}}"
    ))
    label_type <- screen_method
  }

  screen_method <- arg_match(
    screen_method,
    names(ScreenStrategy)
  )
  method_config <- ScreenStrategy[[screen_method]]

  phenotype_class <- arg_match(
    phenotype_class,
    method_config@phenotype_class
  )

  params <- rlang::list2(
    sc_data = sc_data,
    matched_bulk = matched_bulk,
    phenotype = phenotype,
    label_type = label_type,
    phenotype_class = phenotype_class,
    ...
  )
  # * modified params according to specific method
  if (!is.null(method_config@mapper)) {
    params <- method_config@mapper(params)
  }

  res_list <- do.call(method_config@executor, params)
  do.call(ScreenMethodResult, res_list)
}
