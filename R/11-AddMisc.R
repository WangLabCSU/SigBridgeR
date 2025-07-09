# * ------------------- other function ----------------------

#' @title Safely Add Miscellaneous Data to Seurat Object
#'
#' @description
#' Adds arbitrary data to the `@misc` slot of a Seurat object with automatic key
#' conflict resolution. If the key already exists, automatically appends a numeric
#' suffix to ensure unique key naming (e.g., "mykey1", "mykey2").
#'
#' @usage
#' AddMisc(
#'   seurat_obj,
#'   key,
#'   value
#' )
#'
#' @param seurat_obj A Seurat object to modify
#' @param key Character string specifying the name for storing the data. Will be
#'        automatically made unique if already exists in `@misc`.
#' @param value Any R object to store in the `@misc` slot. Common uses include:
#'        - Analysis results (data frames, lists)
#'        - Processing parameters
#'        - Supplemental annotations
#' @param cover Logical indicating whether to overwrite existing data. If
#'        `FALSE` (default).
#'
#' @return The modified Seurat object with added `@misc` data. The original object
#'         structure is preserved with no other modifications.
#'
#' @section Key Generation Rules:
#' 1. If `key` doesn't exist: uses as-is
#' 2. If `key` exists: appends the next available number (e.g., "key1", "key2")
#' 3. If numbered keys exist (e.g., "key2"): increments the highest number
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' seurat_obj <- AddMisc(seurat_obj, "QC_stats", qc_df)
#'
#' # Auto-incrementing example
#' seurat_obj <- AddMisc(seurat_obj, "markers", "markers1")
#' seurat_obj <- AddMisc(seurat_obj, "markers", "markers2")
#' # Stores as "markers1" and "markers2"
#'
#' }
#'
#'
#' @export
#'
AddMisc <- function(seurat_obj, key, value, cover = FALSE) {
  if (key %in% names(seurat_obj@misc) && !cover) {
    existing_keys <- grep(
      paste0("^", key, "\\d*$"),
      names(seurat_obj@misc),
      value = TRUE
    )
    if (length(existing_keys) > 0) {
      nums <- as.integer(sub(paste0("^", key, "(\\d+)$"), "\\1", existing_keys))
      nums <- nums[!is.na(nums)]
      key <- if (length(nums) > 0) {
        paste0(key, max(nums) + 1)
      } else {
        paste0(key, "1")
      }
    } else {
      key <- paste0(key, "1")
    }
  }
  seurat_obj@misc[[key]] <- value

  return(seurat_obj)
}
