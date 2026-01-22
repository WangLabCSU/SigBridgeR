#' @title Add Gene-Level Metadata to Seurat Object (Vectorized, ...-based)
#' @description
#' Add multiple feature-level metadata (vectors or 2D tables) to a Seurat object.
#' Handles vectors (length-checked), matrices, data.frames, data.tables, etc.
#' Columns/variables with duplicated names are suffixed (e.g., "type_1").
#' Gene alignment is auto-detected: rows or cols must match nrow(seurat_obj).
#' @param seurat_obj A Seurat object.
#' @param ... One or more metadata inputs:
#'   - Named/unnamed vectors (length = ngenes)
#'   - 2D objects (matrix, data.frame, data.table, etc.) where
#'     one dimension (rows or cols) has size = ngenes.
#' @param assay Assay name (default: `"RNA"`, fallback to first).
#' @return Modified `seurat_obj` (invisibly).
#'
#' @export
AddMetaFeature <- function(seurat_obj, ..., assay = "RNA") {
  chk::chk_is(seurat_obj, "Seurat")

  assay_names <- names(seurat_obj@assays)
  if (length(assay_names) == 0L) {
    cli::cli_abort(c("x" = "No assays found in Seurat object."))
  }
  if (!(assay %in% assay_names)) {
    assay <- assay_names[[1L]]
  }
  assay_obj <- seurat_obj@assays[[assay]]

  n_genes <- nrow(seurat_obj)
  gene_names <- rownames(seurat_obj)

  # Extract current meta.feature (ensure data.frame)
  mf <- assay_obj[[]]
  if (is.null(mf) || nrow(mf) == 0L) {
    # Create empty data.frame
    mf <- as.data.frame(
      matrix(
        nrow = n_genes,
        ncol = 0L,
        dimnames = list(gene_names, NULL)
      ),
      stringsAsFactors = FALSE
    )
  } else if (!identical(rownames(mf), gene_names)) {
    # Ensure rownames match — critical for later
    cli::cli_abort(c(
      "meta.data features rownames was corrupted and must be reset."
    ))
  }

  dt_mf <- data.table::as.data.table(mf, keep.rownames = "gene")

  # * Process ...

  dots <- rlang::list2(...)
  if (length(dots) == 0L) {
    return(invisible(seurat_obj))
  }

  for (i in seq_along(dots)) {
    meta <- dots[[i]]
    nm <- names(dots)[i]

    # Case 1: Atomic vector (numeric/character/logical/...)
    if (is.atomic(meta) && is.null(dim(meta))) {
      if (length(meta) != n_genes) {
        cli::cli_abort(c(
          "x" = "Input {i} (vector {if (nzchar(nm)) paste0(' \"', nm, '\"')}) has length {.val {length(meta)}}, but object has {.val {n_genes}} genes."
        ))
      }

      colname <- if (nzchar(nm)) nm else paste0("V", i)
      colname <- MakeUniqueName(colname, names(dt_mf))

      data.table::set(dt_mf, j = colname, value = meta)
      next
    }

    # Case 2: 2D object (matrix, data.frame, data.table, DataFrame, etc.)
    dims <- dim(meta)
    if (length(dims) != 2L) {
      cli::cli_abort(c("Input {i} is not a vector or 2D object."))
    }

    nrow_meta <- dims[[1L]]
    ncol_meta <- dims[[2L]]

    # Auto-detect gene dimension: must match n_genes
    gene_dim <- if (nrow_meta == n_genes) {
      1L
    } else if (ncol_meta == n_genes) {
      2L
    } else {
      cli::cli_abort(c(
        "Input {i} has dimensions {.val {nrow_meta}} x {.val {ncol_meta}}; neither matches gene count ({.val {n_genes}})."
      ))
    }

    # Extract as data.table, ensure gene names align
    if (is.matrix(meta)) {
      rn <- rownames(meta) %||% paste0("g", seq_len(nrow_meta))
      cn <- colnames(meta) %||% paste0("V", seq_len(ncol_meta))
      if (gene_dim == 2L) {
        meta <- t(meta)
        rn <- cn
        cn <- rn
      }
      dt_meta <- data.table::as.data.table(meta)
      data.table::setnames(dt_meta, cn)
      data.table::setattr(dt_meta, "row.names", rn)
    } else {
      # data.frame / data.table / tibble
      dt_meta <- data.table::as.data.table(meta, keep.rownames = "gene")

      if (gene_dim == 2L) {
        # genes in columns → melt + dcast is heavy; better transpose via matrix
        # But avoid if possible — prefer user to pass genes in rows
        mat <- data.matrix(dt_meta[, !"gene", with = FALSE])
        rn_orig <- if ("gene" %chin% names(dt_meta)) {
          dt_meta[["gene"]]
        } else {
          rownames(meta)
        }
        cn_orig <- names(dt_meta)[names(dt_meta) != "gene"]

        if (length(rn_orig) != nrow(mat) || length(cn_orig) != ncol(mat)) {
          cli::cli_abort(c(
            "Failed to transpose 2D input {i}: dimension mismatch."
          ))
        }

        mat_t <- Matrix::t(mat)
        dt_meta <- data.table::as.data.table(mat_t)
        data.table::setnames(dt_meta, cn_orig)
        data.table::setattr(dt_meta, "row.names", rn_orig)
      }
    }

    # Now: dt_meta must have rownames = gene names (in some order)
    rn_meta <- rownames(dt_meta)
    if (is.null(rn_meta)) {
      rn_meta <- seq_len(nrow(dt_meta))
    }

    # Reorder to match seurat_obj
    data.table::setkeyv(dt_meta, NULL) # clear keys
    dt_meta <- dt_meta[match(gene_names, rn_meta), ]
    data.table::setattr(dt_meta, "row.names", gene_names)

    # Add each column with unique name
    for (j in seq_along(names(dt_meta))) {
      old_nm <- names(dt_meta)[j]
      new_nm <- MakeUniqueName(old_nm, names(dt_mf))
      if (new_nm != old_nm) {
        cli::cli_warn(
          "Renamed column {.val {old_nm}} to {.val {new_nm}} to avoid duplication."
        )
      }
      # Assign by reference
      data.table::set(dt_mf, j = new_nm, value = dt_meta[[j]])
    }
  }

  # --- Finalize ---------------------------------------------------------------

  # Drop 'gene' column, restore rownames
  rn_final <- dt_mf[["gene"]]
  dt_mf[["gene"]] <- NULL
  mf_final <- as.data.frame(dt_mf, stringsAsFactors = FALSE)
  rownames(mf_final) <- rn_final

  # Assign back
  assay_obj[[]] <- mf_final
  seurat_obj@assays[[assay]] <- assay_obj

  invisible(seurat_obj)
}

#' Helper: make name unique by suffixing _1, _2, ...
#'
#' @keywords internal
MakeUniqueName <- function(name, existing) {
  if (!name %chin% existing) {
    return(name)
  }
  i <- 1L
  candidate <- paste0(name, "_", i)
  while (candidate %chin% existing) {
    i <- i + 1L
    candidate <- paste0(name, "_", i)
  }
  candidate
}
