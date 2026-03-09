#' @title Screen Single-Cell Data Using SCIPAC Algorithm
#' @description
#' Performs single-cell phenotype association analysis using the SCIPAC
#' method. Integrates bulk RNA-seq phenotype information with single-cell transcriptomic data
#' to identify cell populations associated with clinical outcomes.
#'
#' @param matched_bulk Matrix or data frame of preprocessed bulk RNA-seq expression
#'   data (genes × samples). Column names must match names/IDs in \code{phenotype}.
#' @param sc_data A matrix/Matrix (genes × cells) or a Seurat object containing
#'   scRNA-seq data to be screened.
#' @param phenotype Phenotype data, either:
#'   \itemize{
#'     \item Named vector (names match \code{matched_bulk} columns)
#'     \item Data frame with row names matching \code{matched_bulk} columns,
#'           containing \code{"time"} and \code{"status"} columns for survival analysis
#'   }
#' @param label_type Character. Context info specifying the screening method identifier.
#'   Default: \code{"SCIPAC"}.
#' @param phenotype_class Type of phenotypic outcome:
#'   \describe{
#'     \item{\code{"binary"}}{Binary traits (e.g., case/control, responder/non-responder)}
#'     \item{\code{"continuous"}}{Continuous measurements (e.g., age, biomarker levels)}
#'     \item{\code{"survival"}}{Survival information (time-to-event data)}
#'   }
#'   Partial matching is supported.
#' @param hvg Integer. Number of highly variable genes to use for preprocessing.
#'   Default: \code{1000L}.
#' @param do_pca_sc Logical. if \code{TRUE}, first do PCA on \code{sc.dat} and use the rotation matrix on \code{bulk.dat}; if \code{False}, first do PCA on \code{bulk.data} and use the rotation matrix on \code{sc.dat}. The default is \code{FALSE}.
#' @param n_pc Integer. Number of principal components to use. Default: \code{60L}.
#' @param sc_batch_col Character or vector. Batch variable for single-cell data.
#'   If character, should be a column name in \code{sc_data} `@metadata`. If vector,
#'   should be a batch assignment vector matching cell order. Default: \code{NULL}
#'   (no batch correction).
#' @param resolution Integer. Clustering resolution for cell type identification.
#'   Higher values produce more clusters. Default: \code{2L}.
#' @param ela_net_alpha Numeric. Elastic net mixing parameter (0 = ridge, 1 = lasso).
#'   Default: \code{0.4}.
#' @param bt_size Integer. Bootstrap sample size for stability assessment.
#'   Default: \code{50L}.
#' @param ncore Integer. Number of CPU cores for parallel computation.
#'   Default: \code{7L}.
#' @param ci_alpha Numeric. Significance level for confidence intervals.
#'   Default: \code{0.05}.
#' @param nfold Integer. Number of folds for cross-validation for regression models. Default: \code{10L}.
#' @param ... Additional arguments (support `assay`(character), `verbose`(logical) & `seed`(integer)).
#'
#' @return A named list containing:
#'   \describe{
#'     \item{\code{scRNA_data}}{Modified Seurat object with SCIPAC screening results
#'           added as metadata columns. Includes association statistics, p-values,
#'           and confidence intervals for each cell cluster.}
#'     \item{\code{pca_res}}{PCA rotation results from \code{SCIPAC::sc.bulk.pca()}}
#'     \item{\code{cluster_res}}{Clustering results from \code{SCIPAC::seurat.ct()}}
#'   }
#'
#' @section Workflow:
#' \enumerate{
#'   \item \strong{Gene overlap}: Identifies common genes between scRNA-seq and bulk data
#'   \item \strong{Preprocessing}: Selects HVGs and preprocesses both datasets
#'   \item \strong{PCA rotation}: Aligns scRNA-seq and bulk data in shared PCA space
#'   \item \strong{Clustering}: Identifies cell populations using Seurat clustering
#'   \item \strong{Association testing}: Tests cluster-phenotype associations using
#'         elastic net regression with bootstrap validation
#'   \item \strong{Result integration}: Adds screening results to Seurat object metadata
#' }
#'
#' @section Requirements:
#' \itemize{
#'   \item R package: \code{SCIPAC} (install from \code{Exceret/SCIPAC})
#'   \item Gene symbols must be consistent between scRNA-seq and bulk datasets
#'   \item Sufficient overlapping genes (typically >100) for reliable rotation
#' }
#'
#' @section Notes:
#' \itemize{
#'   \item For survival analysis, \code{phenotype} must be a data frame with
#'         \code{"time"} and \code{"status"} columns
#'   \item Binary phenotypes are automatically converted to factors
#'   \item Continuous phenotypes use Gaussian family in elastic net regression
#'   \item Screening parameters are stored in \code{scRNA_data@misc$SCIPAC_params}
#'         for reproducibility
#' }
#'
#' @family screen_method
#' @export
#'
#'
#' @examples
#' \dontrun{
#' # Example 1: Binary phenotype (case/control)
#' library(SCIPAC)
#'
#' # Prepare data
#' bulk_expr <- matrix(rpois(1000, 10), nrow = 200, ncol = 5)
#' colnames(bulk_expr) <- paste0("Sample", 1:5)
#' rownames(bulk_expr) <- paste0("Gene", 1:200)
#'
#' seurat_obj <- CreateSeuratObject(counts = matrix(rpois(2000, 5), nrow = 200, ncol = 10))
#'
#' phenotype_binary <- c(Sample1 = 1, Sample2 = 0, Sample3 = 1, Sample4 = 0, Sample5 = 1)
#'
#' # Run SCIPAC screening
#' result <- DoSCIPAC(
#'   matched_bulk = bulk_expr,
#'   sc_data = seurat_obj,
#'   phenotype = phenotype_binary,
#'   phenotype_class = "binary",
#'   hvg = 100,
#'   resolution = 1,
#'   verbose = TRUE
#' )
#'
#' # Access results
#' screened_seurat <- result$scRNA_data
#' head(screened_seurat[[]])
#'
#' # Example 2: Survival analysis
#' survival_data <- data.frame(
#'   time = c(12.5, 24.3, 18.7, 36.2, 30.1),
#'   status = c(1, 0, 1, 1, 0),
#'   row.names = paste0("Sample", 1:5)
#' )
#'
#' result <- DoSCIPAC(
#'   matched_bulk = bulk_expr,
#'   sc_data = seurat_obj,
#'   phenotype = survival_data,
#'   phenotype_class = "survival",
#'   ela_net_alpha = 0.5,
#'   bt_size = 100,
#'   ncore = 4
#' )
#'
#' # Example 3: Continuous phenotype with custom parameters
#' continuous_pheno <- c(Sample1 = 25.3, Sample2 = 45.7, Sample3 = 38.2,
#'                       Sample4 = 52.1, Sample5 = 31.8)
#'
#' result <- DoSCIPAC(
#'   matched_bulk = bulk_expr,
#'   sc_data = seurat_obj,
#'   phenotype = continuous_pheno,
#'   phenotype_class = "continuous",
#'   hvg = 500,
#'   n_pc = 30,
#'   resolution = 2,
#'   ela_net_alpha = 0.3,
#'   ci_alpha = 0.01,
#'   seed = 42
#' )
#'
#' # View screening parameters used
#' result$scRNA_data@misc$SCIPAC_params
#'
#' # Example 4: With batch correction
#' seurat_obj$batch_id <- sample(c("Batch1", "Batch2"), ncol(seurat_obj), replace = TRUE)
#'
#' result <- DoSCIPAC(
#'   matched_bulk = bulk_expr,
#'   sc_data = seurat_obj,
#'   phenotype = phenotype_binary,
#'   phenotype_class = "binary",
#'   sc_batch_col = "batch_id",
#'   ncore = 8
#' )
#' }
#'
#' @seealso \code{\link[SCIPAC]{SCIPAC}}, \code{\link[SCIPAC]{sc.bulk.pca}},
#'          \code{\link[SCIPAC]{seurat.ct}}, \code{\link{RegisterScreenMethod}}
DoSCIPAC <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = "SCIPAC",
  phenotype_class = c("binary", "survival", "continuous"),
  hvg = 1000L,
  do_pca_sc = FALSE,
  n_pc = 60L,
  sc_batch_col = NULL,
  resolution = 2L,
  ela_net_alpha = 0.4,
  bt_size = 50L,
  ncore = 7L,
  ci_alpha = 0.05,
  nfold = 10L,
  ...
) {
  CheckInstalled("Exceret/SCIPAC")
  # Validate phenotype_class parameter
  purrr::walk(
    .x = c(hvg, n_pc, resolution, bt_size, ncore, nfold),
    .f = chk::chk_integer,
  )
  purrr::walk(
    .x = c(ela_net_alpha, ci_alpha),
    .f = chk::chk_numeric,
  )
  chk::chk_chr(label_type)
  chk::chk_chr(phenotype_class)
  if (!is.null(sc_batch_col)) {
    chk::chk_character(sc_batch_col)
  }
  chk::chk_logical(do_pca_sc)

  family <- switch(
    phenotype_class,
    "binary" = "binomial",
    "survival" = "cox",
    "continuous" = "gaussian",
    cli::cli_abort(c(
      "x" = "phenotype_class must be one of 'binary', 'survival', or 'continuous'"
    ))
  )

  # Extract additional arguments
  dots <- list(...)
  verbose <- dots$verbose %||% getFuncOption("verbose") %||% TRUE
  seed <- dots$seed %||% getFuncOption("seed") %||% 123L
  assay <- dots$assay %||% "RNA"

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green("Start SCIPAC screening"))
    ts_cli$cli_alert_info("Find common variable geme (hvg = {.val {hvg}})")
  }

  sc_bulk_prec <- OverLapSCBulk(
    sc_data = sc_data,
    bulk = matched_bulk,
    hvg = hvg,
    assay = assay
  )
  if (!is.null(sc_batch_col)) {
    sc_batch_col <- if (chk::vld_chr(sc_batch_col)) {
      sc_data[[sc_batch_col]]
    } else {
      sc_batch_col
    }
  }

  if (verbose) {
    ts_cli$cli_alert_info("Running PCA (n_pc = {.val {n_pc}})")
  }

  pca_res <- SCIPAC::sc.bulk.pca(
    sc.dat = sc_bulk_prec$sc_data_preprocessed,
    bulk.dat = sc_bulk_prec$bulk_preprocessed,
    do.pca.sc = do_pca_sc,
    n.pc = n_pc,
    batch_var = sc_batch_col # NULL or a vector
  )

  if (verbose) {
    ts_cli$cli_alert_info("Clustering (resolution = {.val {resolution}})")
  }

  cluster_res <- SCIPAC::seurat.ct(
    sc.dat = pca_res$sc.dat.rot,
    res = resolution
  )

  phenotype <- switch(
    family,
    "binomial" = ,
    "gaussian" = as.factor(phenotype),
    "cox" = as.matrix(phenotype)
  )

  if (verbose) {
    ts_cli$cli_alert_info("Screening ({.val {bt_size} bootstrap})")
  }

  SCIPAC_res <- SCIPAC::SCIPAC(
    bulk.dat = pca_res$bulk.dat.rot,
    y = phenotype,
    family = family,
    ct.res = cluster_res,
    ela.net.alpha = ela_net_alpha,
    bt.size = bt_size,
    numCores = ncore,
    CI.alpha = ci_alpha,
    nfold = nfold
  )

  modified_sc_data <- SeuratObject::AddMetaData(
    object = sc_data,
    metadata = SCIPAC_res
  ) %>%
    AddMisc(
      SCIPAC_para = list(
        phenotype_class = phenotype_class,
        hvg = hvg,
        do_pca_sc = do_pca_sc,
        n_pc = n_pc,
        sc_batch_col = sc_batch_col,
        resolution = resolution,
        ela_net_alpha = ela_net_alpha,
        bt_size = bt_size,
        ncore = ncore,
        ci_alpha = ci_alpha,
        nfold = nfold,
        k = cluster_res$k
      ),
      SCIPAC_type = label_type
    )

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green("SCIPAC screening done"))
  }

  # Return result in expected format
  list(
    scRNA_data = modified_sc_data,
    pca_res = pca_res,
    cluster_res = cluster_res
  )
}


#' @keywords internal
OverLapSCBulk <- function(sc_data, bulk, hvg = 1000L, assay = "RNA") {
  sc_expr <- SeuratObject::LayerData(sc_data, assay = assay, layer = "counts")

  cm_genes <- intersect(rownames(sc_expr), rownames(bulk))

  if (length(cm_genes) == 0L) {
    cli::cli_abort(c(
      "x" = "No overlapping genes between single cell data and bulk data"
    ))
  }

  bulk_new <- as.matrix(bulk[cm_genes, ])

  # * a base matrix
  sc_data_preprocessed <- SCIPAC::obtain.preprocessed.data(
    exprs.data = sc_expr[cm_genes, ],
    hvg = hvg,
    assay = assay
  )
  bulk_preprocessed <- bulk_new[rownames(sc_data_preprocessed), ]

  list(
    "sc_data_preprocessed" = sc_data_preprocessed,
    "bulk_preprocessed" = bulk_preprocessed
  )
}
