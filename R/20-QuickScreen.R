#' @title A Quick Function to Screen Seurat Object
#' @description
#' This function will integrate the matched bulk expression data, phnotype data and the seurat object to filter out cells that are highly correlated with the phenotype, many algorithms are available.
#'
#'
#' @param matched_bulk
#' @param sc_data A Seurat object to be screened.
#' @param phenotype
#' @param label_type A character string indicating the name of phenotype, .
#' @param family This argument is used in the screen_method `Scissor` and `scPAS`, support `gaussian`, `binomial`, `cox`
#' @param phenotype_class This argument is only used in the screen_method `scAB` and `scPP`, supports `binary`, `survival` for `scAB` and supports `binary`, `continuous`, `survival` for `scPP`
#' @param screen_method support `Scissor`, `scPP`, `scPAS`, `scAB`
#' @param ... Other arguments to be passed to the chose screen method.
#'
#' @return A list containing the screened seurat object, in name `scRNA_data`
#' @export
#'
#'
Screen <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type,
  family = c("gaussian", "binomial", "cox"),
  phenotype_class = c("binary", "survival", "continuous"),
  screen_method = c("Scissor", "scPP", "scPAS", "scAB"),
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
          phenotype_class = glue::glue(
            toupper(substr(phenotype_class, 1, 1)),
            tolower(substr(phenotype_class, 2, nchar(phenotype_class)))
          )
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
      "scab" ~
        {
          phenotype_class = tolower(phenotype_class)
          DoscAB(sc_data, matched_bulk, phenotype, label_type, phenotype_class)
        },
      TRUE ~ stop("Screen method not found.")
    )

  return(screened_result)
}
