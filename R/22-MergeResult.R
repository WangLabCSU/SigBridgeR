#' @title Merge Multiple Screening Analysis Results
#' @description
#' Combines results from multiple single-cell screening analyses (e.g., Scissor,
#' scPAS, scPP, scAB) by merging their metadata, slots, and miscellaneous
#' information while preserving the original expression data from the first input.
#' Performs an inner join on cell barcodes to ensure only cells present in all
#' inputs are retained.
#'
#' @param ... Input objects to merge. Can be:
#'   \itemize{
#'     \item Seurat objects
#'     \item Lists containing \code{scRNA_data} (Seurat objects)
#'     \item Mixed combinations of the above
#'   }
#'   The first object serves as the base for merging. When duplicate metadata
#'   columns are found, priority is given to the first occurrence.
#' @param weights Named numeric vector specifying weights for voting when multiple
#'   screening methods produce categorical labels. Names must match screening
#'   method names (case-insensitive). If \code{NULL} (default), no voting is performed.
#' @param ties.method Character. Method for resolving ties in weighted voting:
#'   \describe{
#'     \item{\code{"first"}}{Select the first category alphabetically (default)}
#'     \item{\code{"last"}}{Select the last category alphabetically}
#'     \item{\code{"random"}}{Randomly select among tied categories}
#'   }
#' @param verbose Logical. Whether to print progress messages. Default: inherits
#'   from \code{getOption("SigBridgeR.verbose")}.
#'
#' @return A merged Seurat object containing:
#'   \itemize{
#'     \item Expression data from the first input object
#'     \item Combined metadata from all input objects (inner join on cell IDs)
#'     \item Merged slots (assays, reductions, graphs, images) from all objects
#'     \item Combined miscellaneous information from all objects
#'     \item Optional \code{vote} column in metadata if \code{weight} is provided
#'   }
#'
#' @section Processing Workflow:
#' \enumerate{
#'   \item \strong{Input Validation}: Extracts Seurat objects from inputs (skips invalid objects with warning)
#'   \item \strong{Metadata Extraction}: Converts metadata to data.table format for efficient merging
#'   \item \strong{Cell Intersection}: Retains only cells present in all datasets (inner join)
#'   \item \strong{Weighted Voting}: If \code{weight} provided, performs weighted voting across screening methods
#'   \item \strong{Slot Merging}: Combines assays, reductions, graphs, and images from all objects
#'   \item \strong{Misc Merging}: Aggregates miscellaneous information from all objects
#' }
#'
#' @section Weighted Voting:
#' When \code{weight} is provided, the function:
#' \enumerate{
#'   \item Identifies columns in metadata matching registered screening method names
#'   \item Extracts vote data (must contain "Positive", "Negative", "Neutral", or "Other")
#'   \item Applies \code{\link{WeightedVote}} with specified weights
#'   \item Adds \code{vote} column to merged metadata with final aggregated labels
#' }
#'
#' @section Notes:
#' \itemize{
#'   \item Only cells present in \emph{all} input objects are retained (inner join behavior)
#'   \item Duplicate metadata columns are resolved by keeping the first occurrence
#'   \item For integrating \emph{heterogeneous} single-cell datasets (different samples),
#'         use \code{\link{SCIntegrate}} instead
#'   \item Weights should be non-negative; higher values indicate greater influence in voting
#' }
#'
#' @family result_management
#' @export
#'
#'
#' @examples
#' \dontrun{
#' # Example 1: Merge results from different screening methods
#' scissor_out <- RunScissor(...)
#' scAB_out <- RunscAB(...)
#' scPP_out <- RunscPP(...)
#'
#' merged <- MergeResult(scissor_out, scAB_out, scPP_out)
#'
#' # Check merged metadata
#' head(merged[[]])
#'
#' # Example 2: Weighted voting across methods
#' weights <- c(
#'   Scissor = 3,    # Highest confidence
#'   scAB = 2,       # Medium confidence
#'   scPP = 1        # Lower confidence
#' )
#'
#' merged <- MergeResult(
#'   scissor_out,
#'   scAB_out,
#'   scPP_out,
#'   weight = weights,
#'   ties.method = "first"
#' )
#'
#' # View voting results
#' table(merged$vote)
#'
#' # Example 3: Merge list-containing objects
#' result_list1 <- list(scRNA_data = seurat_obj1, other_info = "...")
#' result_list2 <- list(scRNA_data = seurat_obj2, other_info = "...")
#'
#' merged <- MergeResult(result_list1, result_list2)
#'
#' # Example 4: Check which cells were retained
#' original_cells <- ncol(seurat_obj1)
#' merged_cells <- ncol(merged)
#' cat("Retained", merged_cells, "of", original_cells, "cells\n")
#'
#' # Example 5: Access merged slots
#' # All unique assays from all inputs
#' names(merged@assays)
#'
#' # All unique reductions
#' names(merged@reductions)
#' }
#'
#' @seealso \code{\link{WeightedVote}}, \code{\link{SCIntegrate}}, \code{\link{AddMisc}}
MergeResult <- function(
  ...,
  weights = NULL,
  ties.method = c("random", "first", "last"),
  verbose = SigBridgeRUtils::getFuncOption("verbose")
) {
  args <- rlang::list2(...)
  ..duplicate_cols <- ..vote_cols <- NULL # suppress checking NOTE

  if (length(args) == 0) {
    Abort("Input objects must be provided.")
  }
  # Extract Seurat objects
  seurat_objects <- lapply(args, function(x) {
    if (S7_inherits(x, "ScreenMethodResult")) {
      return(x@scRNA_data)
    } else if (inherits(x, "Seurat")) {
      return(x)
    } else if (is.list(x) && inherits(x$scRNA_data, "Seurat")) {
      return(x$scRNA_data)
    } else {
      cli::cli_warn(
        "Skipping object of class {.code {class(x)}} - not a Seurat object or a list with Seurat object"
      )
      return(NULL)
    }
  })
  seurat_objects <- Filter(Negate(is.null), seurat_objects)

  if (length(seurat_objects) == 0) {
    Abort("No valid Seurat objects found in inputs.")
  }

  # extract metadata
  meta_list <- lapply(seurat_objects, function(x) {
    data.table::as.data.table(x[[]], keep.rownames = "cell_id")
  })

  merged_meta <- Reduce(
    function(x, y) {
      duplicate_cols <- setdiff(intersect(names(x), names(y)), "cell_id")

      if (length(duplicate_cols) > 0) {
        y <- y[, !..duplicate_cols]
      }

      data.table::merge.data.table(x, y, by = "cell_id", all = FALSE)
    },
    meta_list
  )

  # * vote
  if (!is.null(weights)) {
    if (verbose) {
      cli::cli_alert_info("Performing weighted voting")
    }

    existing <- tolower(names(ScreenStrategy))

    method_vote <- tolower(names(weights))
    not_exist <- !method_vote %chin% existing
    if (any(not_exist)) {
      cli::cli_warn(c(
        "x" = "Some voting methods do not exist: {.val {method_vote}}",
        ">" = "These are ignored"
      ))
    }

    vote_cols <- colnames(merged_meta)[
      tolower(colnames(merged_meta)) %chin% method_vote
    ]
    vote_data <- merged_meta[, ..vote_cols]

    weights_needed <- weights[!not_exist]

    vote_label <- WeightedVote(
      vote_data = vote_data,
      weights = weights_needed,
      ties.method = ties.method
    )
    merged_meta$vote <- vote_label
  }

  common_cells <- merged_meta$cell_id
  # check if all cells are present, ignore it if data are identical
  common_cells_len <- length(common_cells)
  first_seurat_cells <- ncol(seurat_objects[[1]])
  if (common_cells_len != first_seurat_cells) {
    cli::cli_warn(
      c(
        "i" = "The {.fun MergeResult} was not originally designed for integrating heterogeneous single-cell datasets",
        ">" = "Only the intersection of cells will be retained",
        ">" = "After intersection, only {.val {common_cells_len}} cells retained from {.val {first_seurat_cells}} cells in the first base object",
        ">" = "Use {.fun SCIntegrate} to integrate heterogeneous single-cell datasets"
      )
    )

    merged_obj <- subset(seurat_objects[[1]], cells = common_cells)
    merged_obj[[]] <- SigBridgeRUtils::Col2Rownames(merged_meta, "cell_id")

    # merge slots
    merged_obj <- Reduce(
      function(merged_obj, slot_type) {
        MergeSlot(
          slot_type = slot_type,
          merged_obj = merged_obj,
          seurat_objects = seurat_objects,
          common_cells = common_cells
        )
      },
      c("assays", "reductions", "graphs", "images"),
      init = merged_obj
    )
  } else {
    # all same, just use the first one as base
    merged_obj <- seurat_objects[[1]]
    merged_obj[[]] <- SigBridgeRUtils::Col2Rownames(merged_meta, "cell_id")
  }

  # merge misc
  all_keys <- unique(unlist(lapply(seurat_objects, function(obj) {
    if (!is.null(obj@misc)) names(obj@misc) else character(0)
  })))

  misc_list <- stats::setNames(vector("list", length(all_keys)), all_keys)

  for (key in all_keys) {
    values <- lapply(seurat_objects, function(obj) {
      if (!is.null(obj@misc) && key %chin% names(obj@misc)) {
        obj@misc[[key]]
      } else {
        NULL
      }
    })

    values <- values[!vapply(X = values, FUN = is.null, FUN.VALUE = logical(1))]

    misc_list[[key]] <- if (length(values) == 1) values[[1]] else values
  }

  merged_obj <- AddMisc(merged_obj, misc_list, cover = TRUE)

  if (verbose) {
    cli::cli_alert_success(
      "Successfully merged {.val {length(seurat_objects)}} objects."
    )
  }

  merged_obj
}

#' @title Helper function to get slot names
#' @description
#' Returns the names of the slots in a Seurat object.
#'
#' @keywords internal
GetSlotNames <- function(obj, slot_type) {
  switch(
    slot_type,
    "assays" = names(obj@assays),
    "graphs" = names(obj@graphs),
    "reductions" = names(obj@reductions),
    "images" = names(obj@images),
    character(0)
  )
}

#' @title Helper function to merge slot
#' @description
#' Merges the slot of a Seurat object (`merged_obj`) with the slot (`slot_type`) of other Seurat objects (`seurat_objects`).
#' `common_cells` is a vector of cell barcodes that are common to all Seurat objects.
#'
#' @keywords internal
MergeSlot <- function(slot_type, merged_obj, seurat_objects, common_cells) {
  base_names <- GetSlotNames(merged_obj, slot_type)
  # The first object is the base object
  for (i in seq_along(seurat_objects)[-1]) {
    current_obj <- seurat_objects[[i]]
    current_names <- GetSlotNames(current_obj, slot_type)

    # Find unique names not in base
    unique_names <- setdiff(current_names, base_names)

    # Add unique items
    for (name in unique_names) {
      item <- switch(
        slot_type,
        "assays" = current_obj@assays[[name]],
        "graphs" = current_obj@graphs[[name]],
        "reductions" = current_obj@reductions[[name]],
        "images" = current_obj@images[[name]]
      )

      # Subset to common cells if applicable
      if (slot_type == "assays") {
        item <- subset(item, cells = common_cells)
        merged_obj@assays[[name]] <- item
      } else if (slot_type == "reductions") {
        valid_cells <- intersect(
          rownames(item),
          common_cells
        )
        if (length(valid_cells) > 0) {
          item <- item[valid_cells, ]
          merged_obj@reductions[[name]] <- item
        }
      } else if (slot_type == "graphs") {
        merged_obj@graphs[[name]] <- item
      } else if (slot_type == "images") {
        merged_obj@images[[name]] <- item
      }
    }

    base_names <- GetSlotNames(merged_obj, slot_type)
  }

  merged_obj
}
