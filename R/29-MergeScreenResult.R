#' @title Merge Multiple Screening Analysis Results
#'
#' @description
#' Combines results from multiple single-cell screening analyses (Scissor, scPAS, scPP, or scAB)
#' by merging their metadata while preserving the original expression data. Performs an inner join
#' on cell barcodes to ensure only cells present in all inputs are retained.
#'
#' @usage
#' MergeResult(...)
#'
#' @param ... Input objects to merge. Can be:
#'        - Seurat objects
#'        - Lists containing `scRNA_data` (Seurat objects)
#'        - Mixed combinations of the above
#'        - The first one will be used as base object for merging, priority given to first one when duplicate columns are found
#'
#' @return A merged Seurat object containing:
#' \itemize{
#'   \item Expression data from the first input object
#'   \item Combined metadata from all input objects
#'   \item Only cells present in all input objects (inner join)
#' }
#'
#' @section Processing Details:
#' 1. Input Validation: Checks for valid Seurat objects or lists containing Seurat objects
#' 2. Metadata Extraction: Collects metadata from all objects
#' 3. Cell Intersection: Retains only cells present in all datasets
#' 4. Object Merging: Creates new Seurat object with combined metadata
#'
#' @examples
#' \dontrun{
#' # Merge two Scissor results
#' merged <- MergeResult(scissor_result1, scissor_result2)
#'
#' # Merge mixed analysis types
#' combined <- MergeResult(scissor_output, scAB_output, scPP_output)
#'
#' # Merge list-containing objects
#' merged_list <- MergeResult(list1, list2, seurat_obj)
#' }
#'
#' @export
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom cli cli_abort cli_alert_warning cli_alert_success
#' @importFrom crayon yellow
#' @importFrom data.table merge.data.table
#' @importFrom purrr map
#'
MergeResult = function(...) {
    args <- list(...)

    if (length(args) == 0) {
        cli::cli_abort(
            "Input objects must be provided."
        )
    }

    seurat_objects <- list()
    seurat_objects <- lapply(args, function(x) {
        if (inherits(x, "Seurat")) {
            return(x)
        } else if (is.list(x) && inherits(x$scRNA_data, "Seurat")) {
            return(x$scRNA_data)
        } else {
            cli::cli_alert_warning(crayon::yellow(
                "Skipping object of class {.code class(x)} - not a Seurat object or list with Seurat object"
            ))
            return(NULL)
        }
    }) %>%
        Filter(Negate(is.null), .)

    if (length(seurat_objects) == 0) {
        cli::cli_abort("No valid Seurat objects found in inputs.")
    }

    meta_list <- lapply(seurat_objects, function(x) {
        data.table::as.data.table(x@meta.data, keep.rownames = "cell_id")
    })

    merged_meta <- Reduce(
        function(x, y) {
            duplicate_cols <- setdiff(intersect(names(x), names(y)), "cell_id")

            if (length(duplicate_cols) > 0) {
                y <- y[, !..duplicate_cols]
            }

            merge(x, y, by = "cell_id", all = FALSE)
        },
        meta_list
    )

    merged_obj <- seurat_objects[[1]]
    merged_obj@meta.data <- merged_meta %>%
        tibble::column_to_rownames("cell_id")

    cli::cli_alert_success(
        "Successfully merged {length(seurat_objects)} objects."
    )

    return(merged_obj)
}
