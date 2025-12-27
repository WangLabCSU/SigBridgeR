# * ---- scRNA-seq preprocessing ----

#' @title Single-Cell RNA-seq Preprocessing Pipeline
#'
#' @description
#' A generic function for standardized preprocessing of single-cell RNA-seq data
#' from multiple sources. Handles data.frame/matrix, AnnData, and Seurat inputs
#' with tumor cell filtering. Implements a complete analysis pipeline
#' from raw data to clustered embeddings.
#'
#' @name SCPreProcess
#'
#' @param sc Input data, one of:
#'    - `data.frame/matrix/dgCMatrix`: Raw count matrix (features x cells)
#'    - `R6`: Python AnnData object, obtained via package `anndata` or `anndataR`
#'    - `Seurat`: Preprocessed Seurat object
#' @param meta_data A data.frame containing metadata for each cell. It will be
#'    added to the Seurat object as `@meta.data`. If sc is an anndata object,
#'    `obs` will be automatically used.
#' @param column2only_tumor A character of column names in `meta_data`, used to
#'    filter the Seurat object to only tumor cells. If `NULL`, no filtering is performed.
#' @param project A character of project name, used to name the Seurat object. Pass to `CreateSeuratObject`.
#' @param min_cells Minimum number of cells that must express a feature for it.
#'    to be included in the analysis. Defaults to `400L`. Pass to `CreateSeuratObject`
#' @param min_features Minimum number of features that must be detected in a
#'    cell for it to be included in the analysis. Defaults to `0L`. Pass to `CreateSeuratObject`
#' @param quality_control A `list` of QC settings. If `NULL`, no QC metrics are computed.
#'    Default: `list(pattern = "^MT-")`.
#'    \describe{
#'      \item{pattern}{A character vector of regex patterns to identify gene groups (e.g., mitochondrial, ribosomal).}
#'    }
#' @param data_filter A `list` of filtering thresholds. If `NULL`, no cell filtering is performed.
#'    Default:
#'    \code{list(
#'      nFeature_RNA_thresh = c(200L, 6000L),
#'      nCount_RNA_thresh = c(500L, 50000L),
#'      percent.mt = 20L,
#'      percent.rp = 60L
#'    )}
#' @param normalization_method Method for normalization: "LogNormalize", "CLR",
#'    or "RC". Defaults to `"LogNormalize"`. Pass to `NormalizeData`
#' @param scale_factor Scaling factor for normalization. Defaults to `10000L`. Pass to `ScaleData`
#' @param scale_features Features to use for scaling. If NULL, uses all variable
#'    features. If `"hvg"`, uses high-variance genes via `VariableFeatures()`. Defaults to `NULL`.
#'    Pass to `ScaleData(features = scale_features)`
#' @param selection_method Method for variable feature selection: "vst", "mvp",
#'    or "disp". Defaults to `"vst"`. Pass to `FindVariableFeatures`
#' @param resolution Resolution parameter for clustering. Higher values lead to
#'    more clusters. Defaults to `0.6`. Pass to `FindClusters`
#' @param dims Dimensions to use for clustering and dimensionality reduction.
#'    If NULL, automatically determined by elbow method. Defaults to `NULL`.
#' @param ... Additional arguments passed to specific methods. Currently supports:
#'    - `verbose`: Logical indicating whether to print progress messages. Defaults to `TRUE`.
#'    - `dims_Neighbors`: Dimensions to use for `FindNeighbors`. Defaults to `NULL`, using `dims`.
#'    - `dims_TSNE`: Dimensions to use for `RunTSNE`. Defaults to `NULL`, using `dims`.
#'    - `dims_UMAP`: Dimensions to use for `RunUMAP`. Defaults to `NULL`, using `dims`.
#'
#' @return A Seurat object containing:
#' \itemize{
#'   \item Normalized and scaled expression data
#'   \item Variable features identified by selection method
#'   \item PCA, t-SNE, and UMAP dimensionality reductions
#'   \item Cluster identities at specified resolution
#'   \item Quality control metrics in `@meta.data`
#'   \item When tumor cells filtered: original dimensions in `@misc$raw_dim`
#'   \item Final dimensions in `@misc$self_dim`
#'   \item Quality control column names in `@misc$qc_colnames`
#' }
#'
#' @details
#' \strong{Quality Control Patterns:}
#' The function supports flexible pattern matching for quality control, for example:
#' \itemize{
#'   \item \code{"^MT-"} - Mitochondrial genes (default)
#'   \item \code{"^RP\[LS\]"} - Ribosomal protein genes
#'   \item \code{"^\[rt\]rna"} - rRNA and tRNA genes
#'   \item Custom patterns using regular expressions
#'   \item Combined patterns: \code{"^MT-|^RP\[LS\]"} for both mitochondrial and ribosomal genes
#' }
#'
#' \strong{Flexible Filtering:}
#' The filtering system dynamically adapts to detected quality control patterns:
#' \itemize{
#'   \item Column names are automatically generated from patterns
#'   \item Multiple thresholds can be specified in \code{data_filter.thresh}
#'   \item Use \code{SigBridgeR:::Pattern2Colname} to determine correct column names for custom patterns if still confused
#' }
#'
#' @examples
#' \dontrun{
#' # Example with matrix input
#' counts_matrix <- matrix(rpois(1000, 5), nrow = 100, ncol = 10)
#' rownames(counts_matrix) <- paste0("Gene", 1:100)
#' colnames(counts_matrix) <- paste0("Cell", 1:10)
#'
#' seurat_obj <- SCPreProcess(
#'   sc = counts_matrix,
#'   project = "TestProject",
#'   min_features = 50,
#'   resolution = 0.8
#' )
#'
#' # Example with tumor cell filtering
#' metadata <- data.frame(
#'   cell_type = c(rep("Tumor", 5), rep("Normal", 5)),
#'   row.names = paste0("Cell", 1:10)
#' )
#'
#' tumor_seurat <- SCPreProcess(
#'   sc = counts_matrix,
#'   meta_data = metadata,
#'   column2only_tumor = "cell_type",
#'   project = "TumorAnalysis"
#' )
#' }
#'
#' @seealso
#' \code{\link[Seurat]{CreateSeuratObject}},
#' \code{\link[Seurat]{NormalizeData}},
#' \code{\link[Seurat]{ScaleData}},
#' \code{\link[Seurat]{FindVariableFeatures}},
#' \code{\link[Seurat]{RunPCA}},
#' \code{\link[Seurat]{RunTSNE}},
#' \code{\link[Seurat]{RunUMAP}},
#' \code{\link[Seurat]{FindNeighbors}},
#' \code{\link[Seurat]{FindClusters}}
#'
NULL


#' @rdname SCPreProcess
#' @export
SCPreProcess <- function(sc, ...) {
  UseMethod("SCPreProcess")
}

#' @rdname SCPreProcess
#' @export
SCPreProcess.default <- function(
  sc,
  meta_data = NULL,
  column2only_tumor = NULL,
  project = "SC_Screening_Proj",
  min_cells = 400L,
  min_features = 0L,
  quality_control = list(
    pattern = c("^MT-")
  ),
  data_filter = list(
    nFeature_RNA_thresh = c(200L, 6000L),
    nCount_RNA_thresh = c(500L, 50000L),
    # * only used when specifed in `quality_control.pattern`
    percent.mt = 20L, # mitochondrial genes
    percent.rp = 60L # ribosomal protein genes
    # ? When combined pattern is used, like `quality_control$pattern <- "^MT-|^RP[LS]"`
    # ? Use `_` to separate different patterns like this:
    # percent.mt_rp = 60L

    # ? When filtering for non-mitochondrial genes and non-ribosomal proteins RNA genes,
    # ? the column names are in lowercase letter form with regular expression symbols removed.
    # `quality_control$pattern <- "^[rt]rna"`
    # Correct threshhold setting is `percent.rt_rna = 60L`

    # ? Use `SigBridgeR:::Pattern2Colname()` to get the correct colname if still confused.
  ),
  normalization_method = "LogNormalize",
  scale_factor = 10000L,
  scale_features = NULL,
  selection_method = "vst",
  resolution = 0.6,
  dims = NULL,
  ...
) {
  if (!is.null(meta_data)) {
    chk::chk_data(meta_data)
  }
  if (!is.null(column2only_tumor)) {
    chk::chk_character(column2only_tumor)
  }

  # dots arguments
  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose") %||% TRUE
  # pass to `FindNeighbors`, `RunTSNE`, `RunUMAP`
  dims_Neighbors <- dots$dims_Neighbors
  dims_TSNE <- dots$dims_TSNE
  dims_UMAP <- dots$dims_UMAP

  sc_seurat <- SeuratObject::CreateSeuratObject(
    counts = sc,
    project = project,
    meta.data = meta_data,
    min.cells = min_cells,
    min.features = min_features
  )

  if (has_pattern(quality_control)) {
    chk::chk_character(quality_control$pattern)
    sc_seurat <- QCPatternDetect(
      obj = sc_seurat,
      pattern = quality_control$pattern,
      verbose = verbose
    )
  }
  if (is_filtering(data_filter)) {
    # first 2 numbers will be used
    chk::chk_lt(
      data_filter$nFeature_RNA_thresh[1],
      data_filter$nFeature_RNA_thresh[2]
    )
    chk::chk_lt(
      data_filter$nCount_RNA_thresh[1],
      data_filter$nCount_RNA_thresh[2]
    )
    sc_seurat <- QCFilter(
      seurat_obj = sc_seurat,
      data_filter.thresh = data_filter,
      verbose = verbose
    )
  }
  sc_seurat <- ProcessSeuratObject(
    obj = sc_seurat,
    normalization_method = normalization_method,
    scale_factor = scale_factor,
    scale_features = scale_features,
    selection_method = selection_method,
    verbose = verbose,
    ...
  )

  sc_seurat <- ClusterAndReduce(
    sc_seurat,
    dims = dims,
    dims_Neighbors = dims_Neighbors,
    dims_TSNE = dims_TSNE,
    dims_UMAP = dims_UMAP,
    resolution = resolution,
    verbose = verbose,
    ...
  )

  FilterTumorCell(
    obj = sc_seurat,
    column2only_tumor = column2only_tumor,
    verbose = verbose
  )
}

# * ---- anndata and anndataR ----

#' @rdname SCPreProcess
#' @export
SCPreProcess.R6 <- function(
  sc,
  meta_data = NULL,
  column2only_tumor = NULL,
  project = "SC_Screening_Proj",
  min_cells = 400L,
  min_features = 0L,
  quality_control = list(
    pattern = c("^MT-")
  ),
  data_filter = list(
    nFeature_RNA_thresh = c(200L, 6000L),
    nCount_RNA_thresh = c(500L, 50000L),
    percent.mt = 20L, # mitochondrial genes
    percent.rp = 60L # ribosomal protein genes
  ),
  normalization_method = "LogNormalize",
  scale_factor = 10000L,
  scale_features = NULL,
  selection_method = "vst",
  resolution = 0.6,
  dims = NULL,
  ...
) {
  if (!is.null(meta_data)) {
    chk::chk_data(meta_data)
  }
  if (!is.null(column2only_tumor)) {
    chk::chk_character(column2only_tumor)
  }

  # dots arguments
  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose") %||% TRUE
  # pass to `FindNeighbors`, `RunTSNE`, `RunUMAP`
  dims_Neighbors <- dots$dims_Neighbors
  dims_TSNE <- dots$dims_TSNE
  dims_UMAP <- dots$dims_UMAP

  # Both `anndata` and `anndataR` are based on R6
  if (is.null(sc$X)) {
    cli::cli_abort(c("x" = "Input must contain $X matrix"))
  }
  if (verbose) {
    cli::cli_text(
      "Start from an anndata object"
    )
  }

  sc_seurat <- if (inherits(sc, "InMemoryAnnData")) {
    # * from AnnDataR
    seurat <- sc$as_Seurat()
    if (!is.null(meta_data)) {
      seurat <- SeuratObject::AddMetaData(seurat, meta_data)
    }
    seurat
  } else if (inherits(sc, "AnnDataR6")) {
    # * from anndata

    meta_data <- if (!is.null(sc$obs) && !is.null(meta_data)) {
      cbind(meta_data, sc$obs)
    } else if (!is.null(sc$obs)) {
      sc$obs
    }
    count_mat <- sc$transpose()$X
    seurat <- SeuratObject::CreateSeuratObject(
      counts = count_mat,
      project = project,
      meta.data = meta_data
    )

    if (ncol(sc$var) != 0) {
      seurat <- AddMetaFeature(seurat, sc$var)
    }
    seurat
  } else {
    cli::cli_abort(c(
      "x" = "Input must be an anndata or anndataR object",
      ">" = "Current input is of class {.cls {class(sc)}}"
    ))
  }
  rm(sc)
  gc(verbose = FALSE)

  if (has_pattern(quality_control)) {
    chk::chk_character(quality_control$pattern)
    sc_seurat <- QCPatternDetect(
      obj = sc_seurat,
      pattern = quality_control$pattern,
      verbose = verbose
    )
  }
  if (is_filtering(data_filter)) {
    # first 2 numbers will be used
    chk::chk_lt(
      data_filter$nFeature_RNA_thresh[1],
      data_filter$nFeature_RNA_thresh[2]
    )
    chk::chk_lt(
      data_filter$nCount_RNA_thresh[1],
      data_filter$nCount_RNA_thresh[2]
    )
    sc_seurat <- QCFilter(
      seurat_obj = sc_seurat,
      data_filter.thresh = data_filter,
      verbose = verbose
    )
  }
  sc_seurat <- ProcessSeuratObject(
    obj = sc_seurat,
    normalization_method = normalization_method,
    scale_factor = scale_factor,
    scale_features = scale_features,
    selection_method = selection_method,
    verbose = verbose,
    ...
  )

  sc_seurat <- ClusterAndReduce(
    sc_seurat,
    dims = dims,
    dims_Neighbors = dims_Neighbors,
    dims_TSNE = dims_TSNE,
    dims_UMAP = dims_UMAP,
    resolution = resolution,
    verbose = verbose,
    ...
  )

  FilterTumorCell(
    obj = sc_seurat,
    column2only_tumor = column2only_tumor,
    verbose = verbose
  )
}


#' @rdname SCPreProcess
#' @export
#'
SCPreProcess.Seurat <- function(
  sc,
  column2only_tumor = NULL,
  ...
) {
  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose")

  if (verbose) {
    ts_cli$cli_alert_info("Start from Seurat object")
  }

  # Validation message can be TRUE or message vector
  valid_msg <- methods::validObject(object = sc, test = TRUE)

  # * Successful validation
  if (is.logical(valid_msg)) {
    return(FilterTumorCell(
      obj = sc,
      column2only_tumor = column2only_tumor,
      verbose = verbose
    ))
  }

  # * Failure to validate the Seurat object
  if (verbose) {
    ts_cli$cli_alert_info(
      "Seurat object validation failed, try repairing it..."
    )
    cli::cli_h3(cli::style_bold("Validation message:"))
    purrr::walk(valid_msg, cli::cli_alert_danger)
    cli::cli_h3(cli::style_bold(
      "Try repairing with {.fn UpdateSeuratObject}:"
    ))
  }

  # Decorator function to wrap the UpdateSeuratObject function
  SafelyUpdateSeuratObject <- purrr::safely(
    SeuratObject::UpdateSeuratObject
  )
  updated_res <- SafelyUpdateSeuratObject(sc)
  # Successful repair
  if (is.null(updated_res$error) && verbose) {
    ts_cli$cli_alert_success(cli::col_green(
      "Successfully repaired the Seurat object"
    ))
  } else if (verbose) {
    # Failure to repair
    cli::cli_alert_danger(updated_res$error)
    cli::cli_warn(
      "Seurat object repair failed. It is recommended to rebuild the Seurat object. Filtering is still being performed but may not be reliable."
    )
  }

  FilterTumorCell(
    obj = sc,
    column2only_tumor = column2only_tumor,
    verbose = verbose
  )
}

# * ---- Preprocess ----

#' @title Process a Seurat object (internal)
#'
#' @description
#' Normalize, find variable features, scale, and run PCA
#'
#' @param obj Seurat object
#' @param normalization_method Normalization method ("LogNormalize", "CLR", or "RC")
#' @param scale_factor Scaling factor for normalization
#' @param scale_features Features to scale
#' @param selection_method Variable feature selection method ("vst", "mvp", or "disp")
#' @param verbose Print progress messages
#' @param ... Other params passed to `Seurat::NormalizeData`, `Seurat::FindVariableFeatures`, and `Seurat::ScaleData`
#' @return Seurat object
#'
#' @keywords internal
#' @family single_cell_preprocess
#'
ProcessSeuratObject <- function(
  obj,
  normalization_method = "LogNormalize",
  scale_factor = 10000,
  scale_features = NULL,
  selection_method = "vst",
  verbose = TRUE,
  ...
) {
  obj <- Seurat::NormalizeData(
    object = obj,
    normalization.method = normalization_method,
    scale.factor = scale_factor,
    verbose = verbose,
    ...
  )
  obj <- Seurat::FindVariableFeatures(
    object = obj,
    selection.method = selection_method,
    verbose = verbose,
    ...
  )

  # 🔍 Resolve scale_features with extended semantics:
  #   • NULL          → scale ALL genes (Seurat default)
  #   • "hvg" / "HVG" → scale VariableFeatures(obj)
  #   • character()   → user-provided gene list
  variable_features <- SeuratObject::VariableFeatures(object = obj, ...)

  obj <- Seurat::ScaleData(
    object = obj,
    features = if (
      is.character(scale_features) &&
        length(scale_features) == 1L &&
        toupper(scale_features) == "HVG"
    ) {
      variable_features
    } else {
      scale_features # NULL
    },
    verbose = verbose,
    ...
  )

  Seurat::RunPCA(
    object = obj,
    features = variable_features,
    verbose = verbose,
    ...
  )
}


#' @title Cluster and reduce dimensionality of single-cell data
#' @description
#' Performs clustering and dimensionality reduction on a single-cell object using PCA components.
#' Automatically determines optimal dimensions when not specified.
#'
#' @param obj A single-cell seurat object containing PCA reductions
#' @param dims Vector of dimensions to use (default: `NULL`, auto-determines using `FindRobustElbow`)
#' @param dims_Neighbors Dimensions for neighbor calculation (default: NULL, uses `dims`)
#' @param dims_TSNE Dimensions for t-SNE  (default: NULL, uses `dims`)
#' @param dims_UMAP Dimensions for UMAP  (default: NULL, uses `dims`)
#' @param resolution Clustering resolution parameter (default: `0.6`)
#' @param verbose Whether to print progress messages (default: `TRUE`)
#' @param ... Additional arguments passed to `FindNeighbors`, `FindClusters`, `RunTSNE`, and `RunUMAP`
#'
#' @return Modified single-cell object with clustering and dimensionality reduction results
#'
#' @keywords internal
#'
#' @note Automatically adjusts input dimensions if they exceed available PCA dimensions
ClusterAndReduce <- function(
  obj,
  run_tsne = TRUE,
  dims = NULL,
  dims_Neighbors = NULL,
  dims_TSNE = NULL,
  dims_UMAP = NULL,
  resolution = 0.6,
  verbose = TRUE,
  ...
) {
  n_pcs <- ncol(obj@reductions$pca)
  if (
    is.null(dims) &&
      any(is.null(dims_Neighbors), is.null(dims_TSNE), is.null(dims_UMAP))
  ) {
    dims <- seq_len(FindRobustElbow(obj, verbose = verbose, ndims = 50))
  }
  if (max(dims) > n_pcs) {
    dims <- seq_len(min(max(dims), n_pcs))
    cli::cli_warn(
      "The input dimension is greater than the dimension in PCA. It is now set to the maximum dimension in PCA."
    )
  }
  dims_Neighbors <- dims_Neighbors %||% dims
  dims_TSNE <- dims_TSNE %||% dims
  dims_UMAP <- dims_UMAP %||% dims

  obj <- Seurat::FindNeighbors(
    object = obj,
    dims = dims_Neighbors,
    verbose = verbose,
    ...
  )
  obj <- Seurat::FindClusters(
    object = obj,
    resolution = resolution,
    verbose = verbose,
    ...
  )

  if (run_tsne) {
    obj <- Seurat::RunTSNE(object = obj, dims = dims_TSNE, ...)
  }

  Seurat::RunUMAP(object = obj, dims = dims_UMAP, verbose = verbose, ...)
}

#' @title Filter tumor cells (internal)
#'
#' @description
#' An internal function that filters tumor cells from a Seurat object based on
#' metadata column values. This function identifies tumor cells using pattern
#' matching on cell type labels and creates a subset containing only tumor cells.
#' It also records dimension information before and after filtering for
#' traceability.
#'
#' @param obj Seurat object with a column to filter out tumor cells.
#' @param column2only_tumor Name of the column to filter out tumor cells.
#' @param verbose logical, whether to print progress messages
#'
#' @return A Seurat object containing only tumor cells, with the following
#'         attributes stored in `@misc`:
#'         \itemize{
#'           \item `self_dim`: Dimensions of the filtered object
#'           \item `raw_dim`: Original dimensions before filtering
#'           \item `column2only_tumor`: The column name used for filtering
#'         }
#'         If `column2only_tumor` is `NULL` or the specified column is not found,
#'         returns the original object unchanged.
#'
#' @keywords internal
#' @family single_cell_preprocess
#'
FilterTumorCell <- function(
  obj,
  column2only_tumor = NULL,
  verbose = TRUE
) {
  raw_dim <- dim(obj)
  obj <- SigBridgeRUtils::AddMisc(obj, self_dim = raw_dim, cover = TRUE)

  if (is.null(column2only_tumor)) {
    return(obj)
  }
  if (!column2only_tumor %chin% colnames(obj[[]])) {
    cli::cli_warn(
      "Column '{.emph {column2only_tumor}}' not found, skip tumor cell filtering"
    )
    return(obj)
  }
  if (verbose) {
    cli::cli_text(
      "Filtering tumor cells with '{.emph {column2only_tumor}}'..."
    )
  }

  labels <- obj[[column2only_tumor]][[1]]
  tumor_cells <- grepl(
    "^[Tt]umo.?r|[Cc]ancer[Mm]alignant|[Nn]eoplasm|[Tt]um",
    labels
  )

  obj <- obj[, tumor_cells]

  SigBridgeRUtils::AddMisc(
    seurat_obj = obj,
    raw_dim = raw_dim,
    self_dim = dim(obj),
    column2only_tumor = column2only_tumor,
    cover = TRUE
  )
}

#' @title Calculate Percentage of Features Matching Patterns
#'
#' @description
#' This function calculates the percentage of counts coming from features matching
#' specified patterns (e.g., mitochondrial genes, ribosomal genes) and adds them
#' as metadata columns to the Seurat object.
#'
#' @param obj A seurat object.
#' @param pattern A character vector or list containing regex patterns to identify mitochondrial
#'    genes, ribosomal protein genes, or other unwanted genes, as well as combinations
#'    of these genes. Customized patterns are supported.
#' @param verbose logical, whether to print progress messages
#' @param ... Additional arguments passed to \code{\link[Seurat]{PercentageFeatureSet}}
#'
#' @details
#' The function automatically generates friendly column names based on the patterns:
#' - "mt" for mitochondrial patterns
#' - "rp" for ribosomal patterns
#' - "rrna" for ribosomal RNA patterns
#' - For combined patterns (using |), creates names like "mt_rp"
#' - For other patterns, creates cleaned lowercase names
#'
#' @family single_cell_preprocess
#' @export
#'
QCPatternDetect <- function(
  obj,
  pattern = c("^MT-", "^mt-", "^RP[SL]", "^MT-|^RP[SL]"),
  verbose = TRUE,
  ...
) {
  if (verbose) {
    cli::cli_text(
      "Using QC patterns {.arg {pattern}} to detect metrics"
    )
  }

  # if `pattern` is a list, convert it to a character vector
  patterns <- unlist(pattern)
  colname_mapping <- stats::setNames(
    paste0("percent.", purrr::map_chr(patterns, Pattern2Colname)),
    patterns
  )
  # Maybe these filter already exist in the object
  existing_cols <- colnames(obj[[]])
  existing <- colname_mapping %chin% existing_cols
  patterns_to_process <- names(colname_mapping)[!existing]

  if (verbose && any(existing)) {
    skipped_cols <- colname_mapping[existing]
    skipped_patterns <- names(skipped_cols)

    purrr::walk2(
      skipped_cols,
      skipped_patterns,
      ~ cli::cli_warn(
        "Column {.val {.x}} already exists. Skipping pattern: {.val {.y}}"
      )
    )
  }

  obj <- purrr::reduce(
    .x = patterns_to_process,
    .f = function(obj_acc, pat) {
      col_name <- colname_mapping[[pat]]
      obj_acc[[col_name]] <- Seurat::PercentageFeatureSet(
        obj_acc,
        pattern = pat,
        ...
      )
      obj_acc
    },
    .init = obj
  )

  # Record these colnames to misc slot for further data filter
  obj@misc$qc_colnames <- unname(colname_mapping)

  obj
}

#' @title Filter Seurat object cells by QC metrics
#'
#' @description
#' Filters cells based on nFeature_RNA and optional QC metrics (e.g. percent.mt, percent.rp),
#' defined in `seurat_obj@misc$qc_colnames` (See \code{\link{QCPatternDetect}}). Only metrics with non-constant, non-all-zero values are used.
#'
#' @param seurat_obj A \code{Seurat} object.
#' @param data_filter.thresh A named list with thresholds. Default:
#'   \code{list(
#'     nFeature_RNA_thresh = c(200L, 6000L),
#'     nCount_RNA_thresh    = c(500L, 25000L),
#'     percent.mt = 20L,
#'     percent.rp = 60L
#'   )}.
#' @param verbose Logical; whether to print progress messages via \code{cli}.
#' @param ... No use
#'
#' @return A filtered \code{Seurat} object.
#' @export
#'
#'
QCFilter <- function(
  seurat_obj,
  data_filter.thresh = list(
    nFeature_RNA_thresh = c(200L, 6000L),
    nCount_RNA_thresh = c(500L, 50000L),
    percent.mt = 20L,
    percent.rp = 60L
  ),
  verbose = TRUE,
  ...
) {
  if (!inherits(seurat_obj, "Seurat")) {
    cli::cli_abort(c("x" = "seurat_obj must be a {.cls Seurat} object."))
  }

  if (!is.list(data_filter.thresh)) {
    cli::cli_abort(c(
      "x" = "{.arg data_filter.thresh} must be a named {.cls list}."
    ))
  }

  if (verbose) {
    cli::cli_text(
      "Filtering cells by {.arg nFeature_RNA}, {.arg nCount_RNA} and {.arg QC metrics}"
    )
  }

  defaults <- list(
    nFeature_RNA_thresh = c(200L, 6000L),
    nCount_RNA_thresh = c(500L, 50000L),
    percent.mt = 20L,
    percent.rp = 60L
  )
  thresh <- utils::modifyList(defaults, data_filter.thresh)

  # Ensure nFeature_RNA_thresh is length-2 integer
  if (length(thresh$nFeature_RNA_thresh) != 2) {
    cli::cli_abort(c(
      "x" = "{.arg nFeature_RNA_thresh} must be a integer vector of length 2."
    ))
  }
  if (length(thresh$nCount_RNA_thresh) != 2) {
    cli::cli_abort(c(
      "x" = "{.arg nCount_RNA_thresh} must be a integer vector of length 2."
    ))
  }
  thresh$nFeature_RNA_thresh <- as.integer(thresh$nFeature_RNA_thresh)
  thresh$nCount_RNA_thresh <- as.integer(thresh$nCount_RNA_thresh)

  # filter expr is a string of the form
  # which is used to subset the Seurat object
  nfeat_condition <- rlang::expr(
    nFeature_RNA > !!thresh$nFeature_RNA_thresh[1] &
      nFeature_RNA < !!thresh$nFeature_RNA_thresh[2]
  )
  ncount_condition <- rlang::expr(
    nCount_RNA > !!thresh$nCount_RNA_thresh[1] &
      nCount_RNA < !!thresh$nCount_RNA_thresh[2]
  )

  # see `QCPatternDetect()` for the column names generation
  qc_colnames <- NULL
  if ("qc_colnames" %in% names(seurat_obj@misc)) {
    qc_colnames <- unlist(seurat_obj@misc$qc_colnames, use.names = FALSE)
  }
  if (is.null(qc_colnames) || length(qc_colnames) == 0) {
    qc_colnames <- character(0)
  }

  # --- Validate QC columns and build conditions -------------------------------
  get_qc_condition <- function(qc_colnames, meta, thresh) {
    if (length(qc_colnames) == 0) {
      return(NULL)
    }

    present_cols <- qc_colnames[qc_colnames %in% colnames(meta)]

    valid_cols <- present_cols[vapply(
      present_cols,
      function(col) {
        x <- meta[[col]]
        any(!is.na(x) & x > 0, na.rm = TRUE) &&
          stats::var(x, na.rm = TRUE) > 0
      },
      logical(1L)
    )]

    if (length(valid_cols) == 0) {
      return(NULL)
    }

    purrr::map(valid_cols, function(col) {
      if (!col %in% names(thresh)) {
        if (verbose) {
          cli::cli_alert_info(
            "No threshold provided for {.val {col}}; skipping."
          )
        }
        return(NULL)
      }
      rlang::expr(!!dplyr::sym(col) < !!thresh[[col]])
    }) %>%
      purrr::compact()
  }
  # be aware of expr is not supported in `subset()`
  meta <- seurat_obj[[]]
  qc_conds <- get_qc_condition(qc_colnames, meta, thresh)

  all_conds <- c(list(nfeat_condition, ncount_condition), qc_conds)
  if (length(all_conds) == 0) {
    cli::cli_abort(c("x" = "No valid filtering conditions generated."))
  }

  full_expr <- purrr::reduce(
    all_conds,
    .f = function(x, y) rlang::expr(!!x & !!y)
  )

  # Evaluate safely
  logical_vec <- base::with(data = meta, expr = rlang::eval_tidy(full_expr))

  if (!is.logical(logical_vec) || length(logical_vec) != nrow(meta)) {
    cli::cli_abort(c(
      "x" = "Internal error: filtering condition did not produce a logical vector of length {.val {nrow(meta)}}."
    ))
  }

  keep_cells <- rownames(meta)[logical_vec]

  if (verbose) {
    n_kept <- length(keep_cells)
    n_total <- nrow(meta)
    pct_off <- if (n_total > 0) 100 * (1 - n_kept / n_total) else 0
    cli::cli_text(sprintf(
      "Kept  %d/%d (%.2f%% off) cells after filtering",
      n_kept,
      n_total,
      pct_off
    ))
  }

  subset(seurat_obj, cells = keep_cells)
}

#' @title convert regex patterns to column names (internal)
#'
#' @description
#' An internal utility function that converts regular expression patterns used
#' for quality control in single-cell RNA-seq analysis into standardized column
#' names for Seurat object metadata. This function handles both single patterns
#' and combined patterns with logical OR operators.
#'
#' @param pat A character string containing a regular expression pattern or
#'            multiple patterns combined with `|` (OR operator).
#'
#' @return A character string with the standardized column name derived from
#'         the input pattern(s). The output follows these conventions:
#'         - Lowercase letters only
#'         - Special characters removed or replaced with underscores
#'         - Common patterns mapped to standardized abbreviations
#'         - Combined patterns sorted alphabetically and joined with underscores
#'
#' @examples
#' \dontrun{
#' # Internal usage examples
#' Pattern2Colname("^MT-")                    # Returns "mt"
#' Pattern2Colname("^RP[LS]")                 # Returns "rp"
#' Pattern2Colname("^[rt]rna")                # Returns "rt_rna"
#' Pattern2Colname("^MT-|^RP[LS]")            # Returns "mt_rp"
#' Pattern2Colname("^HB[AB]?")                # Returns "hb_ab"
#' Pattern2Colname("Custom_Pattern[0-9]+")    # Returns "custom_pattern_0_9"
#' }
#'
#' @keywords internal
#' @family single_cell_preprocess
#'
Pattern2Colname <- function(pat) {
  pat_lower <- tolower(pat)

  if (grepl("\\|", pat_lower)) {
    # Handle combined patterns (with | separator)
    parts <- strsplit(pat_lower, "|", fixed = TRUE)[[1]]
    names <- purrr::map_chr(parts, function(p) {
      dplyr::case_when(
        grepl("mt", p) ~ "mt",
        grepl("rp", p) ~ "rp",
        TRUE ~ tolower(gsub("[^[:alnum:]]", "", p))
      )
    })

    return(paste(sort(unique(names)), collapse = "_"))
  }
  # Handle single patterns
  dplyr::case_when(
    grepl("mt", pat_lower) ~ "mt",
    grepl("rp", pat_lower) ~ "rp",
    TRUE ~ {
      clean <- gsub("[^[:alnum:]]", "_", pat_lower)
      clean <- gsub("_+", "_", clean) # Collapse multiple underscores
      clean <- gsub("^_+|_+$", "", clean) # Trim leading/trailing underscores
      tolower(clean)
    }
  )
}


#' @keywords internal
has_pattern <- function(qc_list) {
  !is.null(qc_list) && !is.null(qc_list$pattern) && length(qc_list$pattern) > 0
}

#' @keywords internal
is_filtering <- function(filter_list) {
  !is.null(filter_list)
}
