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
#' }
#'
#' @importFrom scAB create_scAB select_K scAB findSubset
#' @importFrom cli cli_alert_info
#' @importFrom crayon green
#'
#' @keywords internal
#' @noRd
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
  if (!all(rownames(phenotype) == colnames(matched_bulk))) {
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

  scAB_obj <- create_scAB.v5(
    Object = sc_data,
    bulk_dataset = matched_bulk,
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
    " Run NMF with phenotype and cell-cell similarity regularization at",
    " K = {k}."
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


#' @title scAB_data preprocess
#' @description
#' preprocess the single-cell data, bulk data, and phenotype data.
#'
#' @param Object  Seurat object
#' @param bulk_dataset  matrix of bulk data
#' @param phenotype Phenotype data, a matrix with two columns "time" and "state", or a vector
#' @param method method "survival" or "binary"
#'
#' @return a scAB_data
#' @importFrom preprocessCore normalize.quantiles
#' @export
#'
#' @keywords internal
#' @noRd
#'
create_scAB.v5 <- function(
  Object,
  bulk_dataset,
  phenotype,
  method = c("survival", "binary")
) {
  # cell neighbors
  method = match.arg(method)
  if (is.null(Object@graphs$RNA_snn)) {
    cli::cli_alert_danger(c(
      crayon::red("Error:"),
      "'RNA_snn' is not found, please run FindNeighbors function in Seurat."
    ))
  }
  A <- as.matrix(Object@graphs$RNA_snn)
  diag(A) <- 0
  A[which(A != 0)] <- 1
  degrees <- rowSums(A)
  D <- diag(degrees)
  eps = 2.2204e-256
  D12 <- diag(1 / sqrt(pmax(degrees, eps)))

  L <- D12 %*% (D - A) %*% D12 # Normalized Graph Laplacian
  Dhat <- D12 %*% (D) %*% D12
  Ahat <- D12 %*% (A) %*% D12

  # similarity matrix
  sc_exprs <- as.data.frame(Object@assays$RNA$data)
  common <- intersect(rownames(bulk_dataset), rownames(sc_exprs))
  dataset0 <- cbind(bulk_dataset[common, ], sc_exprs[common, ]) # Dataset before quantile normalization.
  dataset1 <- preprocessCore::normalize.quantiles(as.matrix(dataset0)) # Dataset after  quantile normalization.
  rownames(dataset1) <- rownames(dataset0)
  colnames(dataset1) <- colnames(dataset0)
  Expression_bulk <- dataset1[, 1:ncol(bulk_dataset)]
  Expression_cell <- dataset1[, (ncol(bulk_dataset) + 1):ncol(dataset1)]
  X <- cor(Expression_bulk, Expression_cell)
  X = X / norm(X, "F")

  # phenotype ranking
  if (method == "survival") {
    ss <- scAB::guanrank(phenotype[, c("time", "status")])
    S <- diag(1 - ss[rownames(phenotype), 3])
  } else {
    S <- diag(1 - phenotype)
  }
  # return
  obj <- list(
    X = X,
    S = S,
    L = L,
    D = Dhat,
    A = Ahat,
    phenotype = phenotype,
    method = method
  )
  class(obj) <- "scAB_data"
  return(obj)
}
