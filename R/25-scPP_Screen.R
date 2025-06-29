#
DoscPP = function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = c("Binary", "Continuous", "Survival")
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
  if (length(table(phenotype[2])) == 2 || tolower(label_type) = "binary") {
    gene_list = ScPP::marker_Binary(
      matched_bulk,
      phenotype,
      ref_group = "Normal"
    )
  } else if (
    length(table(phenotype[2])) > 3 || tolower(label_type) == "continuous"
  ) {
    gene_list = ScPP::marker_Continuous(
      matched_bulk,
      phenotype[2],
    )
  } else if (ncol(phenotype) == 3 || tolower(label_type) == "survival") {
    gene_list = ScPP::marker_Survival(
      matched_bulk,
      phenotype,
    )
  } else {
    stop("Unknown phenotype type, please check the `label_type`")
  }

  # *Start screen
  scPP_result <- ScPP::ScPP(sc_data, gene_list)

  sc_meta = scPP_result$metadata %>%
    dplyr::mutate(
      `ScPP` = dplyr::case_when(
        ScPP == "Phenotype+" ~ "Positive",
        ScPP == "Phenotype-" ~ "Negative",
        ScPP == "Background" ~ "Neutral"
      )
    )

  sc_data = sc_data %>% Seurat::AddMetaData(sc_meta$ScPP)

  cli::cli_alert_success(c(
    "[{TimeStamp()}]",
    crayon::green("scPP screening done.")
  ))

  return(sc_data)
}
