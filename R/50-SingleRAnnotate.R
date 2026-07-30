#' @title Annotate Single-Cell Data Using SingleR
#' @family Single_Cell_Annotation_Method
#'
#' @description
#' Performs automated cell type annotation using the \code{SingleR} package. Supports
#' built-in references (e.g., Human Primary Cell Atlas) or user-provided
#' \code{SingleCellExperiment} objects. Integrates seamlessly with \code{Seurat}
#' objects by adding prediction results as metadata columns.
#'
#' @param sc A single-cell dataset, either:
#'   \itemize{
#'     \item A \code{Seurat} object (recommended), or
#'     \item A \code{SingleCellExperiment} object.
#'   }
#' @param verbose Logical. Whether to print progress messages.
#'   Default: inherits from package option \code{getOption("SigBridgeR.verbose")}.
#' @param ensembl Logical scalar indicating whether to convert row names to Ensembl IDs. Genes without a mapping to a non-duplicated Ensembl ID are discarded.
#' @param cell.ont String specifying whether Cell Ontology terms should be included in the colData. If `"nonna"`, all samples without a valid term are discarded; if `"all"`, all samples are returned with (possibly NA) terms; if `"none"`, terms are not added.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub. By default, we use data from the gypsum backend.
#' @param ref Reference data for annotation. One of:
#'   \itemize{
#'     \item \code{"HPCA"}: Uses the \emph{Human Primary Cell Atlas} from the \code{celldex} package.
#'     \item A \code{SingleCellExperiment} object containing pre-labeled reference cells.
#'   }
#'   Note: The placeholder option \code{"custom"} is not valid—users must provide an actual reference object.
#'
#'   A numeric matrix of (usually normalized and log-transformed) expression values from a reference dataset, or a SummarizedExperiment object containing such a matrix; see `SingleR::trainSingleR` for details.
#'   Alternatively, a list or List of SummarizedExperiment objects or numeric matrices containing multiple references. Row names may be different across entries but only the intersection will be used,
#'
#' @param labels A character vector of cell type labels corresponding to the columns of \code{ref}.
#'   Required when \code{ref} is a custom \code{SingleCellExperiment}. For \code{"HPCA"},
#'   defaults to \code{ref$label.main}.
#'
#'   Alternatively, if ref is a list, labels should be a list of the same length. Each element should contain a character vector or factor specifying the labels for the columns of the corresponding element of ref.
#'
#' @param ... Additional arguments passed to \code{SingleR::SingleR()} (e.g., \code{assay.type.ref}, \code{genes}).
#'
#' @return
#'   \itemize{
#'     \item If input is a \code{Seurat} object: returns the same object with three new metadata columns:
#'       \describe{
#'         \item{\code{SingleR_labels}}{Predicted cell type labels.}
#'         \item{\code{SingleR_delta_next}}{Confidence score (difference between top and second-best scores).}
#'         \item{\code{SingleR_pruned_labels}}{Labels after pruning low-confidence assignments.}
#'       }
#'     \item If input is a \code{SingleCellExperiment}: returns the \code{SingleR} prediction result
#'   }
#'
#' @section Requirements:
#'   \itemize{
#'     \item The \code{celldex} package is required when using \code{ref = "HPCA"}.
#'     \item All gene names must be consistent between test and reference datasets (e.g., Ensembl or symbol).
#'   }
#'
#' @examples
#' \dontrun{
#' # Using HPCA reference (requires celldex)
#' annotated_seurat <- SingleRAnnotate(seurat_obj, ref = "HPCA")
#'
#' # Using custom reference
#' library(celldex)
#' ref_sce <- BlueprintEncodeData()
#' annotated_sce <- SingleRAnnotate(
#'   sc = my_sce,
#'   ref = ref_sce,
#'   labels = ref_sce$label.main,
#'   de.n = 30
#' )
#' }
#' @export
SingleRAnnotate <- function(
  sc,
  verbose = getFuncOption("verbose"),
  # * pass to `celldex::HumanPrimaryCellAtlasData`
  ensembl = FALSE,
  cell.ont = c("all", "nonna", "none"),
  legacy = FALSE,
  # * pass to `SingleR::SingleR`
  ref = c("HPCA", "custom"),
  labels = NULL,
  ...
) {
  CheckInstalled("dviraran/SingleR")
  CheckInstalled("celldex", where = "bioc")

  orginally_seurat <- inherits(sc, "Seurat")
  sce_obj <- if (orginally_seurat) {
    Seurat::as.SingleCellExperiment(sc)
  } else {
    sc
  }
  # * Robust check
  chk::chk_length(ref)
  if (ref == "HPCA") {
    # Human Primary Cell Atlas but not the `SingleCellExperiment` object.
    check_installed("celldex")
    ref <- celldex::HumanPrimaryCellAtlasData(
      ensembl = ensembl,
      cell.ont = cell.ont,
      legacy = legacy
    )
    labels <- ref$label.main
  } else if (ref == "custom") {
    # misunderstood the argument `ref`
    Abort("Please specify the reference dataset.")
  } else if (is.null(labels)) {
    # Please find the label from ref data yourself because the `label` slotname may be different from different data.
    # e.g. when set `ref = celldex::HumanPrimaryCellAtlasData()`, there are 3 labels
    Abort("Please specify the `labels` from `ref`.")
  }

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green(
      "[SingleR] Start annotating cell types"
    ))
  }

  # * run
  prediction_sce <- SingleR::SingleR(
    test = sce_obj,
    ref = ref,
    labels = labels,
    ...
  )
  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green("Finish"))
  }

  if (!orginally_seurat) {
    return(prediction_sce)
  }

  sc |>
    Seurat::AddMetaData(prediction_sce$labels, col.name = "SingleR_labels") |>
    Seurat::AddMetaData(
      prediction_sce$delta.next,
      col.name = "SingleR_delta_next"
    ) |>
    Seurat::AddMetaData(
      prediction_sce$pruned.labels,
      col.name = "SingleR_pruned_labels"
    )
}
