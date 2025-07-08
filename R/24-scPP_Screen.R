#' @title Perform scPP screening analysis
#'
#' @description
#' This function performs scPP screening on single-cell data using matched bulk data and phenotype information.
#' It supports binary, continuous, and survival phenotype types.
#'
#' @usage
#' DoscPP(
#'   matched_bulk,
#'   sc_data,
#'   phenotype,
#'   label_type,
#'   phenotype_class = c("Binary", "Continuous", "Survival"),
#'   ref_group = 0,
#'   Log2FC_cutoff = 0.585,
#'   estimate_cutoff = 0.2,
#'   probs = 0.2
#' )
#'
#' @param matched_bulk Bulk expression data (genes × samples) where:
#'        - Column names must match `phenotype` row names
#' @param sc_data Seurat object containing preprocessed single-cell data:
#'        - Normalized counts in `RNA` assay
#' @param phenotype Data frame or tibble or named vector with:
#'        - Rownames matching `matched_bulk` columns
#'        - For survival: must contain time and status columns
#' @param label_type Character specifying phenotype label type (e.g., "SBS1")
#' @param phenotype_class Analysis type (case-sensitive):
#'        - `"Binary"`: Case-control studies (e.g., tumor/normal)
#'        - `"Continuous"`: Quantitative traits (e.g., drug response)
#'        - `"Survival"`: Time-to-event data (requires time/status columns)
#' @param ref_group Reference group for binary comparisons (default: 0)
#' @param Log2FC_cutoff Minimum log2 fold-change for binary markers (default: 0.585)
#' @param estimate_cutoff Effect size threshold for continuous traits (default: 0.2)
#' @param probs Quantile cutoff for cell classification (default: 0.2)
#'
#' @return A list containing:
#' \itemize{
#'   \item{scRNA_data}{Seurat object with added metadata:
#'     \itemize{
#'       \item{ScPP: "Positive"/"Negative"/"Neutral" classification}
#'       \item{label_type: Outcome label used}
#'     }
#'   }
#' }
#'
#' @section Algorithm Steps:
#' 1. Data Validation: Checks sample alignment between bulk and phenotype data
#' 2. Marker Selection: Identifies phenotype-associated genes from bulk data
#' 3. Single-cell Screening: Projects bulk markers onto single-cell data
#' 4. Cell Classification: Categorizes cells based on phenotype association
#'
#' @section Reference:
#' Wang X Lab (2023). "scPP: Single-cell Phenotype Prediction".
#' GitHub: \url{https://github.com/WangX-Lab/ScPP}
#'
#' @examples
#' \dontrun{
#' # Binary phenotype analysis
#' res <- DoscPP(
#'   matched_bulk = bulk_data,
#'   sc_data = seurat_obj,
#'   phenotype = ms_data,
#'   label_type = "SBS1",
#'   phenotype_class = "Binary"
#' )
#'
#' # Survival analysis
#' surv_res <- DoscPP(
#'   sc_data = seurat_obj,
#'   matched_bulk = bulk_data,
#'   phenotype = surv_df,
#'   label_type = "OS_status",
#'   phenotype_class = "Survival"
#' )
#' }
#'
#' @importFrom ScPP marker_Binary marker_Continuous marker_Survival ScPP
#' @import Seurat
#' @import dplyr
#' @importFrom glue glue
#' @importFrom crayon green
#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom tibble rownames_to_column
#'
#' @keywords internal
#' @noRd
#'
DoscPP = function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type,
  phenotype_class = c("Binary", "Continuous", "Survival"),
  ref_group = 1,
  Log2FC_cutoff = 0.585,
  estimate_cutoff = 0.2,
  probs = 0.2
) {
  library(dplyr)
  library(Seurat)

  # robust
  if (!all(rownames(phenotype) == colnames(bulk_dataset))) {
    stop(
      "Please check the rownames of phenotype and colnames of bulk_dataset, they should be the same"
    )
  }

  TimeStamp = function() format(Sys.time(), '%Y/%m/%d %H:%M:%S')

  cli::cli_alert_info(c(
    "[{TimeStamp()}]",
    crayon::green(" Start scPP screening.")
  ))

  cli::cli_alert_info(c(
    "[{TimeStamp()}]",
    " Finding markers..."
  ))

  matched_bulk = as.data.frame(matched_bulk)
  # decide which type of phenotype data is used
  if (is.vector(phenotype)) {
    phenotype = as.data.frame(phenotype) %>%
      tibble::rownames_to_column("Sample") %>%
      dplyr::rename("Feature" := 2)
  }
  if (
    length(table(phenotype[2])) == 2 || tolower(phenotype_class) == "binary"
  ) {
    gene_list = ScPP::marker_Binary(
      bulk_data = matched_bulk,
      features = phenotype,
      ref_group = ref_group,
      Log2FC_cutoff = Log2FC_cutoff
    )
  } else if (
    length(table(phenotype[2])) > 3 || tolower(phenotype_class) == "continuous"
  ) {
    gene_list = ScPP::marker_Continuous(
      bulk_data = matched_bulk,
      features = phenotype[2],
      estimate_cutoff = estimate_cutoff
    )
  } else if (ncol(phenotype) == 3 || tolower(phenotype_class) == "survival") {
    gene_list = ScPP::marker_Survival(
      bulk_data = matched_bulk,
      features = phenotype,
    )
  } else {
    stop(
      "Unknown phenotype type, please check the `phenotype_class` and `phenotype`"
    )
  }

  cli::cli_alert_info(c(
    "[{TimeStamp()}]",
    " Screening..."
  ))

  # *Start screen
  scPP_result <- ScPP::ScPP(sc_data, gene_list, probs = probs)

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
      row.names = rownames(.)
    )

  sc_data <- sc_data %>%
    Seurat::AddMetaData(
      metadata = sc_meta[, c("ScPP", "label_type")],
      col.name = c("ScPP", "label_type")
    )
  
  cli::cli_alert_success(c(
    "[{TimeStamp()}]",
    crayon::green(" scPP screening done.")
  ))

  return(list(scRNA_data = sc_data))
}
