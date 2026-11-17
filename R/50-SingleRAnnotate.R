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
#'   * A `Seurat` object (recommended), or
#'   * A `SingleCellExperiment` object, or
#'   * A [ScreenMethodResult] object, whose `scRNA_data` slot is annotated.
#' @param verbose Logical. Whether to print progress messages.
#'   Default: inherits from package option \code{getOption("SigBridgeR.verbose")}.
#' @param ensembl Logical scalar indicating whether to convert row names to Ensembl IDs. Genes without a mapping to a non-duplicated Ensembl ID are discarded.
#' @param cell.ont String specifying whether Cell Ontology terms should be included in the colData. If `"nonna"`, all samples without a valid term are discarded; if `"all"`, all samples are returned with (possibly NA) terms; if `"none"`, terms are not added.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub. By default, we use data from the gypsum backend.
#' @param ref Reference data for annotation. One of:
#'   * `"HPCA"`: Uses the *Human Primary Cell Atlas* from the `celldex` package.
#'   * A `SingleCellExperiment` object containing pre-labeled reference cells.
#'
#'   Note: The placeholder option `"custom"` is not valid-users must provide an actual reference object.
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
#'
#' * If input is a `Seurat` object: returns the same object with three new metadata columns:
#'   * **`SingleR_labels`**: Predicted cell type labels.
#'   * **`SingleR_delta_next`**: Confidence score (difference between top and second-best scores).
#'   * **`SingleR_pruned_labels`**: Labels after pruning low-confidence assignments.
#' * If input is a `SingleCellExperiment`: returns the `SingleR` prediction result
#' * If input is a [ScreenMethodResult]: returns the annotated `Seurat` object
#'   from its `scRNA_data` slot, as above.
#'
#' @section Requirements:
#'
#' * The `celldex` package is required when using `ref = "HPCA"`.
#' * All gene names must be consistent between test and reference datasets (e.g., Ensembl or symbol).
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
SingleRAnnotate <- function(sc, ...) {
  UseMethod("SingleRAnnotate")
}

#' @rdname SingleRAnnotate
#' @export
SingleRAnnotate.default <- function(sc, ...) {
  supported_classes <- utils::methods("SingleRAnnotate")
  sc_class <- class(sc)[1L]

  Abort(
    "{.fn SingleRAnnotate} doesn't support {.cls {sc_class}} objects.",
    "i" = "Supported classes: {.cls {supported_classes}}.",
    "i" = "Convert your object to a supported class first."
  )
}

#' @rdname SingleRAnnotate
#' @export
`SingleRAnnotate.SigBridgeR::ScreenMethodResult` <- function(sc, ...) {
  SingleRAnnotate.Seurat(sc = sc@scRNA_data, ...)
}

#' @rdname SingleRAnnotate
#' @export
SingleRAnnotate.Seurat <- function(sc, ...) {
  prediction_sce <- SingleRAnnotate.SingleCellExperiment(
    sc = Seurat::as.SingleCellExperiment(sc),
    ...
  )

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

#' @rdname SingleRAnnotate
#' @export
SingleRAnnotate.SingleCellExperiment <- function(
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
  rlang::check_installed("SingleR", action = \(pkg, ...) pak::pak("dviraran/SingleR"))

  # * Robust check
  if (is.character(ref)) {
    ref <- rlang::arg_match(ref)
    if (ref == "custom") {
      # misunderstood the argument `ref`
      Abort(
        "Please specify the reference dataset. Not character `'custom'`",
        type = "[ARG ERROR]"
      )
    }
    # Human Primary Cell Atlas but not the `SingleCellExperiment` object.
    rlang::check_installed("celldex")
    ref <- celldex::HumanPrimaryCellAtlasData(
      ensembl = ensembl,
      cell.ont = rlang::arg_match(cell.ont),
      legacy = legacy
    )
    labels <- ref$label.main
  } else if (is.null(labels)) {
    # Please find the label from ref data yourself because the `label` slotname may be different from different data.
    # e.g. when set `ref = celldex::HumanPrimaryCellAtlasData()`, there are 3 labels
    Abort("Please specify the `labels` from `ref`.", type = "[ARG ERROR]")
  }

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green(
      "[SingleR] Start annotating cell types"
    ))
  }

  # * run
  prediction_sce <- SingleR::SingleR(
    test = sc,
    ref = ref,
    labels = labels,
    ...
  )
  if (verbose) {
    ts_cli$cli_alert_success("Annotation finish")
  }

  prediction_sce
}
