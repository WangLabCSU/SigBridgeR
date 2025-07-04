#' @title Perform scPP screening
#' @description
#' This function performs scPP screening on single-cell data using matched bulk data and phenotype information.
#' It supports binary, continuous, and survival phenotype types.
#'
#' @param matched_bulk A matrix of matched bulk expression data (genes x samples)
#' @param sc_data A Seurat object containing single-cell expression data
#' @param phenotype A data.frame with phenotype information (samples x features)
#' @param phnotype_type Type of phenotype label: "Binary", "Continuous", or "Survival"
#'
#' @return A Seurat object with scPP prediction results added as metadata
#'
#' @details The function first identifies marker genes from bulk data based on phenotype type,
#' then applies scPP screening to single-cell data using these markers. Results are categorized
#' as "Positive", "Negative", or "Neutral" and added to the Seurat object metadata.
#'
#' @examples
#' \dontrun{
#' sc_data <- DoscPP(bulk_data, sc_data, phenotype, "Binary")
#' }
#'
#' @export
DoscPP = function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type,
  phenotype_class = c("Binary", "Continuous", "Survival"),
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
    crayon::green("Start scPP screening.")
  ))

  # decide which type of phenotype data is used
  if (length(table(phenotype[2])) == 2 || tolower(phenotype_class) = "binary") {
    gene_list = ScPP::marker_Binary(
      matched_bulk,
      phenotype,
      ref_group = "Normal"
    )
  } else if (
    length(table(phenotype[2])) > 3 || tolower(phenotype_class) == "continuous"
  ) {
    gene_list = ScPP::marker_Continuous(
      matched_bulk,
      phenotype[2],
    )
  } else if (ncol(phenotype) == 3 || tolower(phenotype_class) == "survival") {
    gene_list = ScPP::marker_Survival(
      matched_bulk,
      phenotype,
    )
  } else {
    stop("Unknown phenotype type, please check the `phenotype_class`")
  }

  # *Start screen
  scPP_result <- ScPP::ScPP(sc_data, gene_list, ...)

  sc_meta = scPP_result$metadata %>%
    dplyr::mutate(
      `ScPP` = dplyr::case_when(
        ScPP == "Phenotype+" ~ glue::glue("Positive"),
        ScPP == "Phenotype-" ~ glue::glue("Negative"),
        ScPP == "Background" ~ glue::glue("Neutral")
      )
    ) %>%
    cbind(
      label_type = label_type,
      row.names = rownames(.) # !needs to test here
    )

  sc_data = sc_data %>% Seurat::AddMetaData(sc_meta$ScPP, sc_meta$label_type)

  cli::cli_alert_success(c(
    "[{TimeStamp()}]",
    crayon::green("scPP screening done.")
  ))

  return(list(scRNA_data = sc_data))
}
