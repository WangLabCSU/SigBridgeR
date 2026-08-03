#' @title Single-Cell Data Screening
#'
#' @md
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
#' - "LP_SGL": see also [DoLP_SGL()]
#' - "PIPET": see also [DoPIPET()]
#' - "SIDISH": see also [DoSIDISH()]
#' - "SCIPAC": see also [DoSCIPAC()]
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
#'   If character, should be a column name in sc_data@metadata. If vector, should
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
#' @section Method Compatibility:
#'
#' | Method | Supported Phenotypes | Additional Parameters |
#' |---|---|---|
#' | Scissor | Binary, Continuous, Survival | alpha, cutoff, path2load_scissor_cache, path2save_scissor_inputs, reliability_test, reliability_test.n, reliability_test.nfold, cell_evaluation, cell_evaluation.benchmark_data, cell_evaluation.FDR, cell_evaluation.bootstrap_n |
#' | scPP | Binary, Continuous, Survival | ref_group, Log2FC_cutoff, estimate_cutoff, probs |
#' | scPAS | Binary, Continuous, Survival | n_components, assay, imputation, nfeature, alpha, network_class, permutation_times, FDR_threshold, independent |
#' | scAB | Binary, Survival | alpha, alpha_2, maxiter, tred |
#' | DEGAS | Binary, Continuous, Survival | sc_data.pheno_colname, select_fraction, tmp_dir, env_params, degas_params, normality_test_method |
#' | LP_SGL | Binary, Continuous, Survival | resolution, alpha, nfold, dge_analysis |
#' | PIPET | Binary, Continuous | group, discretize_method, cutoff, log2FC, p_adjust, show_log2FC, freq_counts, normalize, scale, nPerm, distance |
#' | SIDISH | Survival only | sidish_params, env_params |
#' | SCIPAC | Binary, Continuous, Survival | hvg, do_pca_sc, n_pc, sc_batch_col, resolution, ela_net_alpha, bt_size, ncore, ci_alpha, nfold, assay |
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
    "rSIDISH",
    "SCIPAC"
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
    method_config$phenotype_class
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
  if (!is.null(method_config$mapper)) {
    params <- method_config$mapper(params)
  }

  res_list <- eval(method_config$executor, params)
  exec(ScreenMethodResult, !!!res_list)
}
