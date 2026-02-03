#' @title Annotate Single-Cell Data Using SingleR
#' @keywords Single_Cell_Annotation_Method
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
#'   Default: inherits from package option \code{getOption("SigBridgeRUtils.verbose")}.
#' @param ref Reference data for annotation. One of:
#'   \itemize{
#'     \item \code{"HPCA"}: Uses the \emph{Human Primary Cell Atlas} from the \code{celldex} package.
#'     \item A \code{SingleCellExperiment} object containing pre-labeled reference cells.
#'   }
#'   Note: The placeholder option \code{"Your pre-labelled single cell dataset reference"}
#'   is not valid—users must provide an actual reference object.
#' @param labels A character vector of cell type labels corresponding to the columns of \code{ref}.
#'   Required when \code{ref} is a custom \code{SingleCellExperiment}. For \code{"HPCA"},
#'   defaults to \code{ref$label.main}.
#' @param de.method Method for differential expression during scoring. Passed to \code{SingleR::SingleR()}.
#'   Options include \code{"wilcox"} (default) or \code{"classic"}. See \code{?SingleR::SingleR} for details.
#' @param de.n Integer. Number of top marker genes per label to use in scoring.
#'   Larger values may improve accuracy at the cost of speed. Default: \code{20L}.
#' @param assay.type.test Integer or character. Specifies which assay in \code{sc} to use for testing.
#'   Default: \code{1L} (first assay). Passed to \code{SingleR::SingleR()}.
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
SingleRAnnotate = function(
  sc,
  verbose = getFuncOption("verbose"),
  ref = c("HPCA", "custom"),
  labels = NULL,
  de.method = "wilcox",
  de.n = 20L,
  assay.type.test = 1L,
  ...
) {
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
    rlang::check_installed("celldex")
    ref <- celldex::HumanPrimaryCellAtlasData()
    labels <- ref$label.main
  } else if (ref == "custom") {
    # misunderstood the argument `ref`
    cli::cli_abort(c(
      "x" = "Please specify the reference dataset. It should be a {.cls SingleCellExperiment} object."
    ))
  } else if (is.null(labels)) {
    # Please find the label from ref data yourself because the `label` slotname may be different from different data.
    # e.g. when set `ref = celldex::HumanPrimaryCellAtlasData()`, there are 3 labels
    cli::cli_abort(c(
      "x" = "Please specify the `labels` from `ref`."
    ))
  }
  chk::chk_is(ref, "SingleCellExperiment")
  chk::chk_whole_number(de.n)
  chk::chk_whole_number(assay.type.test)

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green("Start SingleR prediction"))
  }

  # * run
  prediction_sce <- SingleR::SingleR(
    test = sce_obj,
    ref = ref,
    labels = labels,
    de.method = de.method,
    de.n = de.n,
    assay.type.test = assay.type.test,
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
