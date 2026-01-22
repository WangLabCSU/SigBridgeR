#' @title Merge Multiple Screening Analysis Results
#'
#' @description
#' Combines results from multiple single-cell screening analyses (Scissor, scPAS, scPP, or scAB)
#' by merging their metadata and miscellaneous information while preserving the original expression data. Performs an inner join
#' on cell barcodes to ensure only cells present in all inputs are retained.
#'
#' @usage
#' MergeResult(..., verbose = SigBridgeRUtils::getFuncOption("verbose"))
#'
#' @param ... Input objects to merge. Can be:
#'        - Seurat objects
#'        - Lists containing `scRNA_data` (Seurat objects)
#'        - Mixed combinations of the above
#'        - The first one will be used as base object for merging, priority given to first one when duplicate columns are found
#' @param verbose Logical, whether to print a message
#'
#' @return A merged Seurat object containing:
#' \itemize{
#'   \item Expression data from the first input object
#'   \item Combined metadata from all input objects
#'   \item Miscellaneous information from all input objects
#'   \item Only cells present in all input objects (inner join)
#' }
#'
#' @section Processing Details:
#' 1. Input Validation: Checks for valid Seurat objects or lists containing Seurat objects
#' 2. Metadata Extraction: Collects metadata from all objects
#' 3. Cell Intersection: Retains only cells present in all datasets
#' 4. Object Merging: Creates new Seurat object with combined metadata
#' 5. Miscellaneous: Adds miscellaneous information to the merged object
#'
#' @examples
#' \dontrun{
#' # Merge mixed analysis types
#' combined <- MergeResult(scissor_output, scAB_output, scPP_output)
#'
#' # Merge list-containing objects
#' merged_list <- MergeResult(list1, list2, seurat_obj)
#' }
#'
#' @export
#'
MergeResult <- function(
  ...,
  verbose = SigBridgeRUtils::getFuncOption("verbose")
) {
  args <- rlang::list2(...)

  if (length(args) == 0) {
    cli::cli_abort(
      c("x" = "Input objects must be provided."),
      class = "InputsNotFound"
    )
  }
  # Extract Seurat objects
  seurat_objects <- list()
  seurat_objects <- lapply(args, function(x) {
    if (inherits(x, "Seurat")) {
      return(x)
    } else if (is.list(x) && inherits(x$scRNA_data, "Seurat")) {
      return(x$scRNA_data)
    } else {
      cli::cli_warn(
        "Skipping object of class {.code {class(x)}} - not a Seurat object or a list with Seurat object"
      )
      return(NULL)
    }
  }) %>%
    Filter(Negate(is.null), .)

  if (length(seurat_objects) == 0) {
    cli::cli_abort(c("x" = "No valid Seurat objects found in inputs."))
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

  common_cells <- merged_meta$cell_id
  # check if all cells are present, ignore it if data are identical
  common_cells_len <- length(common_cells)
  first_seurat_cells <- ncol(seurat_objects[[1]])
  if (common_cells_len != first_seurat_cells) {
    cli::cli_warn(
      c(
        "i" = "The {.fun MergeResult} was not originally designed for integrating heterogeneous single-cell datasets",
        ">" = " Only the intersection of cells will be retained",
        ">" = "After intersection, only {.val {common_cells_len}} cells retained from {.val {first_seurat_cells}} cells in the first base object"
      )
    )
  }

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

    values <- values[!vapply(values, is.null, FUN.VALUE = logical(1))]

    misc_list[[key]] <- if (length(values) == 1) values[[1]] else values
  }

  merged_obj <- SigBridgeRUtils::AddMisc(merged_obj, misc_list, cover = TRUE)

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
