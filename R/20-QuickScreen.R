#' @title Single-Cell Data Screening
#'
#' @description
#' Integrates matched bulk expression data and phenotype information to identify
#' phenotype-associated cell populations in single-cell RNA-seq data using one of
#' four computational methods. Ensures consistency between bulk and phenotype data
#' before analysis.
#'
#' @param matched_bulk Matrix or data frame of preprocessed bulk RNA-seq expression
#'        data (genes x samples). Column names must match names/IDs in `phenotype`.
#' @param sc_data A Seurat object containing scRNA-seq data to be screened.
#' @param phenotype Phenotype data, either:
#'        - Named vector (names match `matched_bulk` columns), or
#'        - Data frame with row names matching `matched_bulk` columns
#' @param label_type Character specifying phenotype label type (e.g., "SBS1", "time")
#' @param phenotype_class Type of phenotypic outcome (must be consistent with input data):
#'        - `"binary"`: Binary traits (e.g., case/control)
#'        - `"continuous"`: Continuous measurements (only for `Scissor`, `scPAS`, `scPP`)
#'        - `"survival"`: Survival objects
#' @param screen_method Screening algorithm to use, there are four options:
#'        - `"Scissor"`: see also `DoScissor()`
#'        - `"scPP"`: see also `DoscPP()`
#'        - `"scPAS"`: see also `DoscPAS()`
#'        - `"scAB"`: see also `DoscAB()`, no continuous support
#' @param ... Additional method-specific parameters:
#' \describe{
#'   \item{Scissor}{\describe{
#'     \item{scissor_alpha}{(numeric) default 0.05}
#'     \item{scissor_cutoff}{(numeric) default 0.2}
#'     \item{path2load_scissor_cache}{(character) default `NULL`}
#'     \item{dir2save_scissor_inputs}{(character) default `NULl`}
#'     \item{nfold}{(integer) default 10}
#'     \item{reliability_test}{(logical) default FALSE}
#'   }}
#'   \item{scPP}{\describe{
#'     \item{embedding_type}{(character) 嵌入类型，可选"UMAP"或"tSNE"，默认"UMAP"}
#'   }}
#'   \item{scPAS}{\describe{
#'     \item{n_components}{(integer) 主成分数量，默认10}
#'   }}
#'   \item{scAB}{\describe{
#'     \item{bandwidth}{(numeric) 核密度估计带宽，默认0.1}
#'   }}
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{scRNA_data}{Filtered Seurat object with phenotype-associated cells}
#'   \item{matched_samples}{Vector of samples used in the analysis}
#'   \item{method_output}{Method-specific output objects}
#' }
#'
#'
#' @section Data Matching Requirements:
#' - `matched_bulk` column names and `phenotype` names/rownames must be identical
#' - Phenotype values must correspond to bulk samples (not directly to single cells)
#' - Mismatches will trigger an error before analysis begins
#'
#' @section Method Compatibility:
#'
#' | **Method** | **Supported Phenotypes**      | **Additional Parameters**      |
#' |------------|-------------------------------|---------------------------------|
#' | `Scissor`  | All three types               | `alpha`, `lambda`               |
#' | `scPP`     | All three types               | `embedding_type`                |
#' | `scPAS`    | All three types               | `n_components`                  |
#' | `scAB`     | Binary/Survival               | `bandwidth`                     |
#'
#'
#'
#' @seealso Associated functions:
#' \itemize{
#'   \item \code{\link{DoScissor}}
#'   \item \code{\link{DoscPP}}
#'   \item \code{\link{DoscPAS}}
#'   \item \code{\link{DoscAB}}
#' }
#'
#'
#' @export
#' @import dplyr
#' @importFrom glue glue
#'
Screen <- function(
    matched_bulk,
    sc_data,
    phenotype,
    label_type,
    phenotype_class = c("binary", "survival", "continuous"),
    screen_method = c("Scissor", "scPP", "scPAS", "scAB"),
    ...
) {
    library(dplyr)
    if (length(screen_method) != 1) {
        cli::cli_abort(c("x" = "Only one {.arg screen_method} is allowed."))
    }

    phenotype_class = tolower(phenotype_class)
    if (length(phenotype_class) != 1) {
        cli::cli_abort(c("x" = "Only one {.arg phenotype_class} is allowed."))
    } else if (!phenotype_class %in% c("binary", "survival", "continuous")) {
        cli::cli_abort(c(
            "x" = "Invalid {.arg phenotype_class = {phenotype_class}}.",
            "i" = " Must be one of {.val binary}, {.val survival}, or {.val continuous}."
        ))
    }

    screened_result = tolower(screen_method) %>%
        switch(
            "scissor" = {
                family = switch(
                    phenotype_class,
                    "binary" = "binomial",
                    "survival" = "cox",
                    "continuous" = "gaussian"
                )

                DoScissor(
                    sc_data = sc_data,
                    matched_bulk = matched_bulk,
                    phenotype = phenotype,
                    label_type = label_type,
                    scissor_family = family, # "gaussian", "binomial", "cox"
                    ...
                )
            },
            "scpas" = {
                family = switch(
                    phenotype_class,
                    "binary" = "binomial",
                    "survival" = "cox",
                    "continuous" = "gaussian",
                )

                DoscPAS(
                    sc_data = sc_data,
                    matched_bulk = matched_bulk,
                    phenotype = phenotype,
                    label_type = label_type,
                    scPAS_family = family, # "gaussian", "binomial", "cox"
                    ...
                )
            },
            "scpp" = {
                phenotype_class = glue::glue(
                    toupper(substr(phenotype_class, 1, 1)),
                    tolower(substr(phenotype_class, 2, nchar(phenotype_class)))
                )

                DoscPP(
                    sc_data = sc_data,
                    matched_bulk = matched_bulk,
                    phenotype = phenotype,
                    label_type = label_type,
                    phenotype_class = phenotype_class, # "Binary", "Continuous", "Survival"
                    ...
                )
            },
            "scab" = {
                if (phenotype_class == "continuous") {
                    stop("scAB does not support continuous phenotype.")
                }

                DoscAB(
                    sc_data = sc_data,
                    matched_bulk = matched_bulk,
                    phenotype = phenotype,
                    label_type = label_type,
                    phenotype_class = phenotype_class, # "Binary", "Survival"
                    ...
                )
            },
            TRUE ~ stop("Screen method not found.")
        )

    return(screened_result)
}
