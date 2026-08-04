#' @title Add Gene-Level Metadata to Seurat Object (Vectorized, ...-based)
#' @description
#' Add multiple feature-level metadata (vectors or 2D tables) to a Seurat object.
#' Handles vectors (length-checked), matrices, data.frames, data.tables, etc.
#' Columns/variables with duplicated names are suffixed (e.g., "type_1").
#' Gene alignment is auto-detected: rows or cols must match nrow(seurat_obj).
#' @param seurat_obj A Seurat object.
#' @param ... One or more metadata inputs:
#'   - Named/unnamed vectors (length = ngenes)
#'   - 2D objects (matrix, data.frame, etc.) where
#'     one dimension (rows or cols) has size = ngenes.
#' @param assay Assay name (default: `"RNA"`, fallback to first).
#' @return Modified `seurat_obj` (invisibly).
#'
#' @export
AddMetaFeature <- function(seurat_obj, ..., assay = "RNA") {
  chk::chk_is(seurat_obj, "Seurat")

  assay_names <- names(seurat_obj@assays)
  if (!(assay %chin% assay_names)) {
    assay <- SeuratObject::DefaultAssay(seurat_obj)
  }
  assay_obj <- Seurat::GetAssay(object = seurat_obj, assay = assay)
  feature_names <- rownames(assay_obj)

  dots <- list2(...)
  if (length(dots) == 0L) {
    cli::cli_warn("No inputs were provided.")
    return(seurat_obj)
  }
  dots_names <- names2(dots)

  for (i in seq_along(dots)) {
    col_name_i <- dots_names[[i]]
    metadata_i <- dots[[i]]

    if (!nzchar(col_name_i) && is.vector(dots[[i]])) {
      Abort(
        "Column name is empty, please provide a name for the metadata.",
        "E.g, col_name = rep('type', {length(feature_names)})"
      )
    }
    if (is.null(names(metadata_i))) {
      cli::cli_warn(
        "Metadata is not named, trying to use feature names to fill the names."
      )
      if (length(feature_names) != length(metadata_i)) {
        Abort(
          "Metadata length does not match feature length, please check the input"
        )
      }
      names(metadata_i) <- feature_names
    }

    assay_obj <- SeuratObject::AddMetaData(
      object = assay_obj,
      metadata = metadata_i,
      col.name = if (is_2d(metadata_i)) NULL else col_name_i
    )
  }
  seurat_obj[[assay]] <- assay_obj
  invisible(seurat_obj)
}
