#' @title Single-Cell Data Screening
#'
#' @description
#' Integrates matched bulk expression data and phenotype information to identify
#' phenotype-associated cell populations in single-cell RNA-seq data using one of
#' four computational methods. Ensures consistency between bulk and phenotype data
#' before analysis.
#'
#' @param matched_bulk Matrix or data frame of preprocessed bulk RNA-seq expression
#'        data (genes x samples). Column names must match names/IDs in `phenotype`.
#' @param sc_data A matrix/Matrix (genes x cells) or a Seurat object containing scRNA-seq data to be screened.
#' @param phenotype Phenotype data, either:
#'        - Named vector (names match `matched_bulk` columns), or
#'        - Patient survival Data frame with row names matching `matched_bulk` columns, colnames named "time" and "status"
#' @param label_type Character specifying phenotype label type (e.g., "SBS1", "time")
#' @param phenotype_class Type of phenotypic outcome (must be consistent with input data):
#'        - `"binary"`: Binary traits (e.g., case/control)
#'        - `"continuous"`: Continuous measurements
#'        - `"survival"`: Survival infomation
#' @param screen_method Screening algorithm to use, there are seven options:
#' - "Scissor": see also DoScissor()
#' - "scPP": see also DoscPP()
#' - "scPAS": see also DoscPAS()
#' - "scAB": see also DoscAB(), continuous phenotype is not supported
#' - "DEGAS": see also DoDEGAS()
#' - "LP_SGL": see also DoLP_SGL()
#' - "PIPET": see also DoPIPET(), using survival as phenotype is not supported
#' @param ... Additional method-specific parameters:
#' \describe{
#' \item{Scissor}{\describe{
#' \item{alpha}{(numeric or NULL) Significance threshold. When NULL, alpha will keep increasing iteratively until the corresponding cells are screened out, default 0.05}
#' \item{cutoff}{(numeric) A threshold for terminating the iteration of alpha, only work when alpha is NULL, default 0.2}
#' \item{path2load_scissor_cache}{(character) default NULL}
#' \item{path2save_scissor_inputs}{(character) A path to save the intermediary data. By using path2load_scissor_cache, the intermediary data can be loaded from the specified path. default "Scissor_inputs.RData"}
#' \item{reliability_test}{(logical) Whether to perform reliability test, default FALSE}
#' \item{reliability_test.nfold}{(integer) Cross-validation folds for reliability test, default 10}
#' \item{reliability_test.n}{(integer) Number of cells to use for reliability test, default 10}
#' \item{cell_evaluation}{(logical) Whether to perform cell evaluation, default FALSE}
#' \item{cell_evaluation.benchmark_data}{.RData Benchmark data for cell evaluation, default NULL}
#' \item{cell_evaluation.FDR}{(numeric) FDR threshold for cell evaluation, default 0.05}
#' \item{cell_evaluation.bootstrap_n}{(integer) Number of bootstrap samples for cell evaluation, default 10}
#' }}
#' \item{scPP}{\describe{
#' \item{ref_group}{(integer or character) Reference group or baseline for binary comparisons, e.g. "Normal" for Tumor/Normal studies and 0 for 0/1 case-control studies. default: 0}
#' \item{Log2FC_cutoff}{(numeric) Minimum log2 fold-change for binary markers, default 0.585}
#' \item{estimate_cutoff}{(numeric) Effect size threshold for continuous traits, default 0.2}
#' \item{probs}{(numeric) Quantile cutoff for cell classification, default 0.2}
#' }}
#' \item{scPAS}{\describe{
#' \item{assay}{(character) Assay to use from sc_data, default "RNA"}
#' \item{imputation}{(logical) Whether to perform imputation, default FALSE}
#' \item{nfeature}{(integer) Number of features to select, default 3000}
#' \item{alpha}{(numeric or NULL) Significance threshold, When NULL, alpha will keep increasing iteratively until the corresponding cells are screened out, default 0.01}
#' \item{independent}{(logical) The background distribution of risk scores is constructed independently of each cell. default: TRUE}
#' \item{network_class}{(character) Network class to use. default: 'SC', indicating gene-gene similarity networks derived from single-cell data. The other one is 'bulk'.}
#' \item{permutation_times}{(integer) Number of permutations, default 2000}
#' \item{FDR_threshold}{(numeric) FDR value threshold for identifying phenotype-associated cells default 0.05}
#' }}
#' \item{scAB}{\describe{
#' \item{alpha}{(numeric) Coefficient of phenotype regularization ,default 0.005}
#' \item{alpha_2}{(numeric) Coefficent of cell-cell similarity regularization, default 0.005}
#' \item{maxiter}{(integer) NMF optimization iterations, default 2000}
#' \item{tred}{(integer) Z-score threshold, default 2}
#' }}
#' \item{DEGAS}{\describe{
#' \item{sc_data.pheno_colname}{(character) Phenotype column name in sc_data, default "NULL"}
#' \item{select_fraction}{(numeric) Fraction of cells to select for DEGAS, default 0.05}
#' \item{tmp_dir}{(character) Temporary directory for DEGAS, default "NULL"}
#' \item{env_params}{(list) Environment parameters for DEGAS, default "list()"}
#' \item{degas_params}{(list) DEGAS parameters, default "list()"}
#' \item{normality_test_method}{(character) Normality test method for DEGAS, default "jarque-bera"}
#' }}
#' \item{LP_SGL}{\describe{
#' \item{resolution}{(numeric) Resolution parameter for Leiden clustering, default 0.6}
#' \item{alpha}{(numeric) Alpha parameter for SGL balancing L1 and L2 penalties, default 0.5}
#' \item{nfold}{(integer) Number of folds for cross-validation, default 5}
#' \item{dge_analysis}{(list) Differential expression analysis settings:
#' \itemize{
#' \item{run: (logical) Whether to run DEG analysis, default FALSE}
#' \item{logFC_threshold: (numeric) Log fold change threshold, default 1}
#' \item{pval_threshold: (numeric) P-value threshold, default 0.05}
#' }}
#' \item{PIPET}{\describe{
#' \item{group}{(character or NULL) Name of a metadata column (e.g., `"orig.ident"`) to stratify cells before screening. When `NULL` (default), screening is performed globally across all cells.}
#' \item{discretize_method}{(character) Strategy to binarize continuous phenotypes internally before marker identification. One of:
#'   \itemize{
#'     \item{\code{"median"} (default):} Equivalent to 2-quantile split (i.e., median threshold).
#'     \item{\code{"kmeans"}:} Two-cluster k-means on the continuous phenotype.
#'     \item{\code{"custom"}:} User-defined cutoffs via `cutoff`.
#'   }}
#' \item{cutoff}{(numeric vector or NULL) Required only if \code{discretize_method = "custom"}. Specifies interior breakpoints on the *normalized, log2-transformed phenotype scale* (i.e., after `scale(log2(x + 1))`). Must be sorted ascending and of length `n_group - 1`.}
#' \item{label_type}{(character) Phenotype label type (e.g., `"PIPET_SBS1"`), stored in `scRNA_data@misc`. Default: `"PIPET"`.}
#' \item{log2FC}{(numeric) Absolute log2 fold-change cutoff for differential expression marker selection in bulk data (via DESeq2-like analysis). Default: 1.}
#' \item{p_adjust}{(numeric) Adjusted p-value (FDR) cutoff for marker gene selection. Default: 0.05.}
#' \item{show_log2FC}{(logical) Whether to annotate markers with signed log2FC direction (e.g., `CD3D_up`). Default: \code{TRUE}.}
#' \item{freq_counts}{(integer or NULL) Minimum number of cells a gene must be expressed in to be retained in scRNA-seq data preprocessing. Default: \code{NULL} (no filtering).}
#' \item{normalize}{(logical) Whether to apply log-normalization (`LogNormalize`) to scRNA-seq counts prior to correlation. Default: \code{TRUE}.}
#' \item{scale}{(logical) Whether to scale (center + unit-variance) gene expression across cells before computing distances. Default: \code{TRUE}.}
#' \item{nPerm}{(integer) Number of label permutations to assess significance of correlation scores. Default: 1000.}
#' \item{distance}{(character) Distance or similarity metric for template matching. Supported:
#'   \code{"cosine"} (default),
#'   \code{"pearson"},
#'   \code{"spearman"},
#'   \code{"kendall"},
#'   \code{"euclidean"},
#'   \code{"maximum"}.}
#' \item{seed}{(integer or NULL) Random seed for reproducibility in marker creation and permutation tests. Default: inherits from `getFuncOption("seed")`.}
#' \item{verbose}{(logical) Whether to print progress messages. Default: inherits from `getFuncOption("verbose")`.}
#' \item{parallel}{(logical) Whether to enable parallel permutations (requires `future::plan()` pre-set). Default: \code{FALSE}.}
#' }}
#' }
#' }}
#'
#' @return A list containing:
#' \describe{
#' \item{scRNA_data}{A Seurat object with phenotype-associated cells labelled in `meta.data` column}
#' \item{Some screen_result}{Important information about the screened result related to the selected method}
#' }
#'
#'
#' @section Data Matching Requirements:
#' - matched_bulk column names and phenotype names/rownames must be identical
#' - Phenotype values must correspond to bulk samples (not directly to single cells)
#' - Mismatches will trigger an error before analysis begins, and there is a built-in pre-run check.
#'
#' @section Method Compatibility:
#'
#' | Method | Supported Phenotypes | Additional Parameters |
#' |------------|-------------------------------|---------------------------------|
#' | Scissor | All three types | alpha, cutoff, path2load_scissor_cache, path2save_scissor_inputs, reliability_test, reliability_test.n,reliability_test.nfold, cell_evaluation,cell_evaluation.benchmark_data,cell_evaluation.FDR,cell_evaluation.bootstrap_n |
#' | scPP | All three types | ref_group, Log2FC_cutoff, estimate_cutoff, probs |
#' | scPAS | All three types | n_components ,assay, imputation,nfeature, alpha,network_class,permutation_times,FDR_threshold,independent |
#' | scAB | Binary/Survival | alpha, alpha_2, maxiter, tred |
#' | DEGAS | All three types | sc_data.pheno_colname,select_fraction,tmp_dir,env_params,degas_params,normality_test_method |
#' | LP_SGL | All three types | resolution, alpha, nfold, dge_analysis |
#' | PIPET | Binary/Continuous | group, discretize_method, cutoff, log2FC, p_adjust, show_log2FC, freq_counts, normalize, scale, nPerm, distance |
#'
#'
#' @seealso Associated functions:
#' \itemize{
#' \item \code{\link{DoScissor}}
#' \item \code{\link{DoscPP}}
#' \item \code{\link{DoscPAS}}
#' \item \code{\link{DoscAB}}
#' \item \code{\link{DoDEGAS}}
#' \item \code{\link{DoLP_SGL}}
#' \item \code{\link{DoPIPET}}
#' }
#'
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
    "PIPET"
  ),
  ...
) {
  chk::chk_length(phenotype_class)
  chk::chk_length(screen_method)

  if (is.null(label_type) || length(label_type) != 1) {
    cli::cli_alert_info(c(
      "i" = "{.var label_type} not specified or not of length 1, using {.val {screen_method}}"
    ))
    label_type <- screen_method
  }
  available_phenotype_class <- c("binary", "survival", "continuous")
  available_screen_method <- c(
    "Scissor",
    "scPP",
    "scPAS",
    "scAB",
    "DEGAS",
    "LP_SGL",
    "PIPET"
  )
  phenotype_class <- SigBridgeRUtils::MatchArg(
    phenotype_class,
    available_phenotype_class,
    NULL
  )
  screen_method <- SigBridgeRUtils::MatchArg(
    screen_method,
    available_screen_method,
    NULL
  )

  switch(
    screen_method,
    "Scissor" = {
      family <- switch(
        phenotype_class,
        "binary" = "binomial",
        "survival" = "cox",
        "continuous" = "gaussian"
      )

      DoScissor(
        sc_data = sc_data,
        matched_bulk = matched_bulk,
        phenotype = phenotype,
        label_type = label_type,
        scissor_family = family, # "gaussian", "binomial", "cox"
        ...
      )
    },
    "scPAS" = {
      family <- switch(
        phenotype_class,
        "binary" = "binomial",
        "survival" = "cox",
        "continuous" = "gaussian",
      )

      DoscPAS(
        sc_data = sc_data,
        matched_bulk = matched_bulk,
        phenotype = phenotype,
        label_type = label_type,
        scPAS_family = family, # "gaussian", "binomial", "cox"
        ...
      )
    },
    "scPP" = {
      phenotype_class <- glue::glue(
        toupper(substr(phenotype_class, 1, 1)),
        tolower(substr(phenotype_class, 2, nchar(phenotype_class)))
      )

      DoscPP(
        sc_data = sc_data,
        matched_bulk = matched_bulk,
        phenotype = phenotype,
        label_type = label_type,
        phenotype_class = phenotype_class, # "Binary", "Continuous", "Survival"
        ...
      )
    },
    "scAB" = {
      if (phenotype_class == "continuous") {
        cli::cli_abort(c(
          "x" = "[{.fun Screen}]: {.strong scAB} does not support continuous phenotype."
        ))
      }

      DoscAB(
        sc_data = sc_data,
        matched_bulk = matched_bulk,
        phenotype = phenotype,
        label_type = label_type,
        phenotype_class = phenotype_class, # "Binary", "Survival"
        ...
      )
    },
    "DEGAS" = {
      DoDEGAS(
        sc_data = sc_data,
        matched_bulk = matched_bulk,
        phenotype = phenotype,
        label_type = label_type,
        phenotype_class = phenotype_class, # "Binary", "Survival", "Continuous"
        ...
      )
    },
    "LP_SGL" = {
      LPSGL_family = switch(
        phenotype_class,
        "binary" = 'logit',
        "survival" = 'cox',
        "continuous" = 'linear'
      )
      DoLP_SGL(
        sc_data = sc_data,
        matched_bulk = matched_bulk,
        phenotype = phenotype,
        label_type = label_type,
        LPSGL_family = LPSGL_family, # "Binary", "Survival", "Continuous"
        ...
      )
    },
    "PIPET" = {
      if (phenotype_class == "survival") {
        cli::cli_abort(c(
          "x" = "[{.fun Screen}]: {.strong PIPET} does not support survival phenotype."
        ))
      }
      DoPIPET(
        matched_bulk = matched_bulk,
        sc_data = sc_data,
        phenotype = phenotype,
        phenotype_class = phenotype_class,
        ...
      )
    }
  )
}
