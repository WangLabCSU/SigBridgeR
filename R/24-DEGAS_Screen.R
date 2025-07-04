# !----ATTENTION: A model must be trained as benchmark----

#
DoDEGAS = function(matched_bulk, sc_data, phenotype, label_type) {
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

  return(list(scRNA_data = 1))
}
