#' @title A Quick Function to Screen Seurat Object
#' @description
#' This function will integrate the matched bulk expression data, phnotype data and the seurat object to filter out cells that are highly correlated with the phenotype, many algorithms are available.
#'
#'
#' @param matched_bulk
#' @param sc_data
#' @param phenotype
#' @param label_type
#' @param family
#' @param phenotype_class
#' @param screen_method
#' @param ...
#'
#' @return
#' @export
#'
#'
Screen <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type,
  family = c("gaussian", "binomial", "cox"),
  phenotype_class = c("Binary", "Continuous", "Survival"),
  screen_method = c("Scissor", "scPP", "scPAS", "DEGAS", "scAB"),
  ...
) {
  library(dplyr)
  if (length(screen_method) > 1) {
    stop("Only one screen method is allowed.")
  }

  screened_seurat = tolower(screen_method) %>%
    dplyr::case_when(
      "scissor" ~
        {
          family = tolower(family)
          if (
            length(family) != 1 && !family %in% c("gaussian", "binomial", "cox")
          ) {
            stop("`family` must be one of gaussian, binomial, cox")
          }

          DoScissor(
            sc_data = sc_data,
            matched_bulk = matched_bulk,
            phenotype = phenotype,
            label_type = label_type,
            scissor_family = family # "gaussian", "binomial", "cox"
          )
        },
      "scpas" ~
        {
          family = tolower(family)
          if (
            length(family) != 1 && !family %in% c("gaussian", "binomial", "cox")
          ) {
            stop("`family` must be one of gaussian, binomial, cox")
          }

          DoscPAS(
            sc_data = sc_data,
            matched_bulk = matched_bulk,
            phenotype = phenotype,
            label_type = label_type,
            scPAS_family = family # "gaussian", "binomial", "cox"
          )
        },
      "scpp" ~
        {
          if (
            length(phenotype_class) != 1 &&
              !phenotype_class %in% c("Binary", "Continuous", "Survival")
          ) {
            stop(
              "`phenotype_class` must be one of Binary, Continuous, Survival"
            )
          }
          DoscPP(
            sc_data = sc_data,
            matched_bulk = matched_bulk,
            phenotype = phenotype,
            label_type = label_type,
            phenotype_class = phenotype_class # "Binary", "Continuous", "Survival"
          )
        },
      "degas" ~ DoDEGAS(sc_data, matched_bulk, phenotype, label_type, family),
      "scab" ~ DoscAB(sc_data, matched_bulk, phenotype, label_type, family),
      TRUE ~ stop("Screen method not found.")
    )

  return(screened_result)
}
