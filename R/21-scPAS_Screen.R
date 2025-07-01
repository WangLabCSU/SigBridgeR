# ---- 2. Do scPAS ----

#' @title Perform scPAS Screening Analysis
#' @description
#' This function performs scPAS screening analysis by integrating bulk and single-cell RNA-seq data.
#' It includes data filtering steps and wraps the core scPAS::scPAS function.
#'
#' @param matched_bulk Bulk RNA-seq data (genes x samples)
#' @param sc_data Single-cell RNA-seq data (Seurat object and preprocessed)
#' @param phenotype Phenotype data frame with sample annotations
#' @param label_type Column name in phenotype containing sample labels
#' @param assay Assay to use from sc_data (default: 'RNA')
#' @param imputation Logical, whether to perform imputation (default: FALSE)
#' @param nfeature Number of features to select (default: 3000, indicating that the top 3000 highly variable genes are selected for model training
#' @param alpha Significance threshold (default: 0.01)
#' @param gene_RNAcount_filter Minimum gene expression threshold (default: 20)
#' @param bulk_0_filter_thresh Maximum proportion of zeros allowed in bulk data (default: 0.25)
#' @param network_class Network class to use (default: 'SC', indicating gene-gene similarity networks derived from single-cell data.)
#' @param family Model family for analysis (options: "cox", "gaussian", "binomial")
#'
#' @return A Seurat object from scPAS analysis
#' @export
DoscPAS = function(
  matched_bulk,
  sc_data,
  phenotype,
  ms_type,
  assay = 'RNA',
  imputation = F,
  nfeature = 3000,
  alpha = 0.01,
  extra_filter = FALSE,
  gene_RNAcount_filter = 20,
  bulk_0_filter_thresh = 0.25,
  network_class = 'SC',
  family = c("cox", "gaussian", "binomial"),
  ...
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
    crayon::green("Start scPAS screening.")
  ))

  if (extra_filter) {
    # *sc filter
    keep_genes <- rownames(sc_data)[
      Matrix::rowSums(sc_data@assays$RNA@counts > 1) >= gene_RNA_count_filter
    ]
    sc_data <- subset(sc_data, features = keep_genes)
    # *bulk filter zero data
    matched_bulk <- matched_bulk %>%
      Matrix::as.matrix() %>%
      .[
        apply(., 1, function(x) sum(x == 0) < bulk_0_filter_thresh * ncol(.)),
      ]
  } else {
    matched_bulk = as.matrix(matched_bulk)
  }

  scPAS_result <- scPAS::scPAS(
    bulk_dataset = matched_bulk,
    sc_dataset = sc_data,
    assay = 'RNA',
    tag = ms_type,
    phenotype = phenotype,
    imputation = imputation,
    nfeature = nfeature,
    alpha = alpha,
    network_class = network_class,
    family = family,
    ...
  )

  # *rename level
  scPAS_result@meta.data = scPAS_result@meta.data %>%
    dplyr::mutate(
      `scPAS` = dplyr::case_when(
        scPAS == "scPAS+" ~ "Positive",
        scPAS == "scPAS-" ~ "Negative",
        scPAS == "0" ~ "Neutral"
      )
    ) %>%
    cbind(
      metadata = data.frame(ms_type = ms_type),
      row.names = colnames(.)
    )

  cli::cli_alert_success(c(
    "[{TimeStamp()}]",
    crayon::green("scPAS screening done.")
  ))

  return(scPAS_result)
}
