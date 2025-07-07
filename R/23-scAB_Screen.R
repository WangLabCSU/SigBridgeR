#' @title Perform scAB Screening Analysis
#'
#' @description
#' Implements the scAB algorithm to identify phenotype-associated cell subpopulations
#' in single-cell RNA-seq data by integrating matched bulk expression and phenotype
#' information. Uses non-negative matrix factorization (NMF) with dual regularization
#' for phenotype association and cell-cell similarity.
#'
#' @param matched_bulk Normalized bulk expression matrix (genes × samples) where:
#'        - Columns match `phenotype` row names
#'        - Genes match features in `sc_data`
#' @param sc_data Seurat object containing preprocessed single-cell data:
#' @param phenotype Data frame with clinical annotations where:
#'        - Rows correspond to `matched_bulk` columns
#'        - For survival: contains `time` and `status` columns
#' @param label_type Character specifying phenotype label type (e.g., "SBS1", "time")
#' @param phenotype_class Analysis mode:
#'        - `"binary"`: Case-control design (e.g., responder/non-responder)
#'        - `"survival"`: Time-to-event analysis data.frame
#' @param alpha Coefficient of phenotype regularization (default=0.005).
#' @param alpha_2 Coefficent of cell-cell similarity regularization (default=5e-05).
#' @param maxiter NMF optimization iterations (default=2000).
#' @param tred Z-score threshold (default=2).
#'
#' @return A list containing:
#' \describe{
#'   \item{scRNA_data}{Filtered Seurat object with selected cells}
#'   \item{scAB_result}{scAB screening result}
#'
#' @section Reference:
#' Zhang Q, Jin S, Zou X (2022). "scAB detects multiresolution cell states with
#' clinical significance by integrating single-cell genomics and bulk sequencing data."
#'
#' @examples
#' \dontrun{
#' # Binary phenotype example
#' result <- DoscAB(
#'   matched_bulk = bulk_matrix,
#'   sc_data = seurat_obj,
#'   phenotype = clinical_df,
#'   label_type = "disease_status",
#'   phenotype_class = "binary",
#'   alpha = 0.005,
#'   alpha_2 = 5e-05,
#'   maxiter = 2000,
#'   tred = 2
#' )
#'
#' # Survival phenotype example
#' surv_result <- DoscAB(
#'   matched_bulk = bulk_matrix,
#'   sc_data = seurat_obj,
#'   phenotype = survival_df,
#'   label_type = "OS_status",
#'   phenotype_class = "survival",
#'   maxiter = 3000
#' )
#' }
#'
#' @export
#' @importFrom scAB create_scAB select_K scAB findSubset
#' @importFrom cli cli_alert_info
#' @importFrom crayon green
#'
DoscAB <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type,
  phenotype_class = c("binary", "survival"),
  alpha = 0.005,
  alpha_2 = 5e-05,
  maxiter = 2000,
  tred = 2
) {
  library(dplyr)

  # robust
  if (!all(rownames(phenotype) == colnames(bulk_dataset))) {
    stop(
      "Please check the rownames of phenotype and colnames of bulk_dataset, they should be the same"
    )
  }
  TimeStamp = function() format(Sys.time(), '%Y/%m/%d %H:%M:%S')

  cli::cli_alert_info(c(
    "[{TimeStamp()}]",
    crayon::green(" Start scAB screening.")
  ))

  sc_data = sc_data %>%
    Seurat::AddMetaData(rep(label_type, ncol(.)), col.name = "label_type")

  # `Obejct` cannot be changed to `Object`
  scAB_obj <- scAB::create_scAB(
    Obejct = sc_data,
    sc_data = matched_bulk,
    phenotype = phenotype,
    method = phenotype_class
  )

  cli::cli_alert_info(c(
    "[{TimeStamp()}]",
    " Selecting K"
  ))

  k <- scAB::select_K(scAB_obj)

  cli::cli_alert_info(c(
    "[{TimeStamp()}]",
    "Run NMF with phenotype and cell-cell similarity regularization"
  ))

  scAB_result <- scAB::scAB(
    Object = scAB_obj,
    K = k,
    alpha = alpha,
    alpha_2 = alpha_2,
    maxiter = maxiter
  )

  cli::cli_alert_info(c(
    "[{TimeStamp()}]",
    " Screening cells"
  ))

  sc_data <- scAB::findSubset(sc_data, scAB_Object = scAB_result, tred = tred)

  sc_data@meta.data <- sc_data@meta.data %>%
    dplyr::rename(scAB = scAB_select) %>%
    dplyr::mutate(
      scAB = case_when(
        scAB == "Other cells" ~ "Other",
        scAB == "scAB+ cells" ~ "Positive",
        TRUE ~ NA
      )
    )

  cli::cli_alert_info(c(
    "[{TimeStamp()}]",
    crayon::green(" scAB screening done.")
  ))

  return(list(scRNA_data = sc_data, scAB_result = scAB_result))
}
