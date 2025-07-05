#' @title scAB Screening
#' @description
#'
#'
#' @param matched_bulk
#' @param sc_data
#' @param phenotype
#' @param label_type
#' @param phenotype_class
#' @param tred
#' @param alpha
#' @param alpha_2
#' @param maxiter
#' @param tred
#'
#' @return
#' @export
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

  cli::cli_alert_info(c(
    "[{TimeStamp()}]",
    crayon::green(" scAB screening done.")
  ))

  return(list(scRNA_data = sc_data, scAB_result = scAB_result))
}
