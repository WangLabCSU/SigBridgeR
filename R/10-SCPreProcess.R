# * ---- scRNA-seq preprocessing ----

#' @title Single-Cell RNA-seq Preprocessing Pipeline
#' @name SCPreProcess
#'
#' @description
#' A unified interface for standardized single-cell RNA-seq preprocessing.
#' It accepts raw counts (matrix/data.frame), AnnData (via `anndata` or `anndataR`), or existing Seurat objects.
#' The analysis flow is fully customizable via a character string \code{pipeline} and a configuration list \code{params}.
#'
#' @param sc Input data. Can be:
#' \itemize{
#'   \item \strong{Matrix/Data frame}: Raw count matrix (genes x cells).
#'   \item \strong{AnnData}: Python AnnData object (read via \code{anndata} or \code{anndataR} packages).
#'   \item \strong{Seurat}: A Seurat object (automatically validated and repaired if necessary).
#' }
#' @param pipeline A character string defining the processing steps and order.
#' Characters map to Seurat functions:
#' \itemize{
#'   \item \code{'o'}: \code{CreateSeuratObject} (Must be the first step and cannot be deleted)
#'   \item \code{'n'}: \code{NormalizeData}
#'   \item \code{'s'}: \code{ScaleData}
#'   \item \code{'v'}: \code{FindVariableFeatures}
#'   \item \code{'p'}: \code{RunPCA}
#'   \item \code{'e'}: \code{FindNeighbors} (Because "n" is used)
#'   \item \code{'c'}: \code{FindClusters}
#'   \item \code{'t'}: \code{RunTSNE}
#'   \item \code{'u'}: \code{RunUMAP}
#'   \item \code{'r'}: \code{SCTransform} (Alternative to n/s/v)
#' }
#' Default is \code{"onsvpcetu"}.
#' @param params A named list of lists containing arguments for each pipeline step.
#' Keys match the pipeline characters (e.g., \code{params$n} for \code{NormalizeData}).
#' Default structure:
#' \preformatted{
#' list(
#'   o = list(project = "SC_Screen_Proj", min.cells = 400L), # do not pass `counts`
#'   n = list(),             # NormalizeData args
#'   s = list(),             # ScaleData args
#'   v = list(),             # FindVariableFeatures args
#'   c = list(resolution = 0.6),
#'   ...
#' )
#' }
#' @param quality_control A list containing regex patterns for QC metric calculation. (See [QCPatternDetect])
#' Default: \code{list(pattern = "^MT-")}. Detected metrics (e.g., percent.mt) are added to meta.data.
#' @param data_filter A list of thresholds for cell filtering.
#' Default: \code{nFeature_RNA} (200-6000), \code{nCount_RNA} (500-50000), \code{percent.mt} (<20), \code{percent.rp} (<60).
#' Only metrics detected via \code{quality_control} are filtered, i.e., nFeature_RNA, nCount_RNA and percent.mt.
#' @param column2only_tumor Optional character. Column name in metadata to filter for tumor cells
#' (matches "Tumor", "Cancer", "Malignant", etc.). If \code{NULL}, no filtering is performed.
#' @param ... Additional arguments for backward compatibility (mapped to \code{params}) or verbose control.
#'
#' @return A processed Seurat object with reductions, clusters, and QC metrics.
#'
#' @details
#' \strong{Pipeline Strategy:}
#' The function parses the \code{pipeline} string and executes corresponding Seurat functions in order.
#' To use \code{SCTransform}, simply change the pipeline string (e.g., \code{"orpetu"}) and provide parameters in \code{params$r}.
#'
#' \strong{Quality Control & Filtering:}
#' QC metrics are generated based on regex patterns in \code{quality_control}.
#' Cells are then filtered based on thresholds in \code{data_filter}.
#' Column names for filtering are auto-generated (e.g., pattern "^MT-" -> filter "percent.mt").
#' If confused about the column name, use \code{SigBridgeR:::Pattern2Colname()}.
#'
#' @examples
#' \dontrun{
#' # 1. Standard pipeline (LogNormalize -> Scale -> PCA -> UMAP)
#' obj <- SCPreProcess(
#'   sc = counts_matrix,
#'   pipeline = "onsvpcetu",
#'   params = list(c = list(resolution = 0.8))
#' )
#'
#' # 2. SCTransform pipeline
#' obj_sct <- SCPreProcess(
#'   sc = counts_matrix,
#'   pipeline = "orpcu", # Create -> SCT -> PCA -> Clusters -> UMAP
#'   quality_control = list(pattern = c("^MT-", "^RP[LS]")),
#'   params = list(
#'     r = list(vars.to.regress = "percent.mt")
#'   )
#' )
#'
#' # 3. Start from AnnData with tumor filtering
#' adata_object <- anndataR::read_h5ad("data.h5ad")
#' obj_ad <- SCPreProcess(
#'   sc = adata_object,
#'   column2only_tumor = "tissue"
#' )
#' }
#'
#' @family single_cell_preprocess
#' @export
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
  ...,
  pipeline = "onsvpcetu",
  params = list(
    # * CreateSeuratObject
    o = list(
      project = "SC_Screen_Proj",
      min.cells = 400L
    ),
    # * NormalizeData
    n = list(),
    # * ScaleData
    s = list(),
    # * FindVariableFeatures
    v = list(),
    # * RunPCA
    p = list(),
    # * FindNeightbors
    e = list(),
    # * FindClusters
    c = list(
      resolution = 0.6
    ),
    # * RunTSNE
    t = list(),
    # * RunUMAP
    u = list()
    # * SCTransform
    # r = list()
  ),
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
  column2only_tumor = NULL
) {
  # * check
  if (!is.null(column2only_tumor)) {
    chk::chk_character(column2only_tumor)
  }
  chk::chk_character(pipeline)
  purrr::walk(c(quality_control, data_filter, params), ~ chk::chk_list)
  purrr::walk(params, ~ chk::chk_list)

  # * dots arguments, to be compatible last version
  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose") %||% TRUE

  params <- compatible_with_3.0.2(..., params = params) # backward compatibility
  params <- handle_usr_params(params) # combine with default params

  # * handling pipeline
  steps <- unlist(strsplit(pipeline, ""))
  steps_to_run <- steps[steps %chin% names(sc_processing_strategies)] # a vector

  unknown <- setdiff(steps, steps_to_run)
  if (length(unknown) != 0) {
    cli::cli_warn(
      "[{.fun SCIntegrate}]: Unknown pipeline steps: {.val {unknown}}"
    )
  }
  if (steps_to_run[[1]] != "o") {
    cli::cli_abort(
      "[{.fun SCPreProcess}]: The first step of {.arg pipeline} must be 'o' for CreateSeuratObject"
    )
  }
  if (!is.null(params$o$counts)) {
    cli::cli_abort(c(
      "x" = "[{.fun SCPreProcess}]: The parameter {.arg params$o$counts} is deprecated, please use {.arg sc} instead"
    ))
  }

  sc_seurat <- rlang::exec(
    SeuratObject::CreateSeuratObject,
    counts = sc,
    !!!params$o
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

  # a list containing functions
  execution_queue <- sc_processing_strategies[steps_to_run]

  # the first step is CreateSeuratObject("o")
  for (step in 2:length(steps_to_run)) {
    # step: an index
    step_fun <- execution_queue[[step]]
    sc_seurat <- step_fun(
      object = sc_seurat,
      # steps_to_run[[step]] -- a letter
      params = params[[steps_to_run[[step]]]]
    )
  }

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
  ...,
  pipeline = "onsvpcetu",
  params = list(
    # * CreateSeuratObject
    o = list(
      project = "SC_Screen_Proj",
      min.cells = 400L
    ),
    # * NormalizeData
    n = list(),
    # * ScaleData
    s = list(),
    # * FindVariableFeatures
    v = list(),
    # * RunPCA
    p = list(),
    # * FindNeightbors
    e = list(),
    # * FindClusters
    c = list(
      resolution = 0.6
    ),
    # * RunTSNE
    t = list(),
    # * RunUMAP
    u = list()
    # * SCTransform
    # r = list()
  ),
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
  column2only_tumor = NULL
) {
  # Both `anndata` and `anndataR` are based on R6

  # * check
  if (!is.null(column2only_tumor)) {
    chk::chk_character(column2only_tumor)
  }
  chk::chk_character(pipeline)
  purrr::walk(c(quality_control, data_filter, params), ~ chk::chk_list)
  purrr::walk(params, ~ chk::chk_list)

  # * dots arguments, to be compatible last version
  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose") %||% TRUE

  params <- compatible_with_3.0.2(..., params = params) # backward compatibility
  params <- handle_usr_params(params) # combine with default params

  # * handling pipeline
  steps <- unlist(strsplit(pipeline, ""))
  steps_to_run <- steps[steps %chin% names(sc_processing_strategies)] # a vector

  unknown <- setdiff(steps, steps_to_run)
  if (length(unknown) != 0) {
    cli::cli_warn(
      "[{.fun SCPreProcess}]: Unknown pipeline steps: {.val {unknown}}"
    )
  }
  if (steps_to_run[[1]] != "o") {
    cli::cli_abort(
      "[{.fun SCPreProcess}]: The first step of {.arg pipeline} must be 'o' for CreateSeuratObject"
    )
  }
  if (!is.null(params$o$counts)) {
    cli::cli_abort(c(
      "x" = "[{.fun SCPreProcess}]: The parameter {.arg params$o$counts} is deprecated, please use {.arg sc} instead"
    ))
  }

  if (is.null(sc$X)) {
    cli::cli_abort(c(
      "x" = "[{.fun SCPreProcess}]: {.arg sc} must contain $X matrix"
    ))
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
    if (params$o$min.cells != 0) {
      gene_cell_counts <- SigBridgeRUtils::rowSums3(
        SeuratObject::LayerData(seurat, layer = "counts") > 0
      )
      seurat <- seurat[gene_cell_counts >= min_cells, ]
    }
    if (params$o$min.features != 0) {
      cell_gene_counts <- SigBridgeRUtils::colSums3(
        SeuratObject::LayerData(seurat, layer = "counts") > 0
      )
      seurat <- seurat[, cell_gene_counts >= min_features]
    }
    seurat
  } else if (inherits(sc, "AnnDataR6")) {
    # * from anndata

    params$o$meta.data <- if (!is.null(sc$obs)) {
      dplyr::bind_cols(params$o$meta.data, sc$obs)
    } else if (!is.null(sc$obs)) {
      params$o$meta.data
    }
    count_mat <- sc$transpose()$X

    if (any(count_mat != floor(count_mat))) {
      cli::cli_warn(
        "[{.fun SCPreProcess}]: Input count matrix in {.arg sc} contains {.emph non-integer} values"
      )
    }

    seurat <- rlang::exec(
      SeuratObject::CreateSeuratObject,
      counts = count_mat,
      !!!params$o
    )

    if (ncol(sc$var) != 0) {
      seurat <- AddMetaFeature(seurat, sc$var)
    }
    seurat
  } else {
    cli::cli_abort(c(
      "x" = "{.arg sc} must be an anndata or anndataR object",
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

  # a list containing functions
  execution_queue <- sc_processing_strategies[steps_to_run]

  # the first stop is CreateSeuratObject("o")
  # step: an index
  for (step in 2:length(steps_to_run)) {
    # step_fun: a function
    step_fun <- execution_queue[[step]]
    sc_seurat <- step_fun(
      object = sc_seurat,
      # steps_to_run[[step]]: a letter
      params = params[[steps_to_run[[step]]]]
    )
  }

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
    cli::cli_text("Start from Seurat object")
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
    cli::cli_alert_info("Seurat object validation failed, try repairing it...")
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
    cli::cli_alert_success(cli::col_green(
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

#' @keywords internal
GetVars2Regress <- function(seurat_obj, verbose = TRUE) {
  if (!"qc_colnames" %in% names(seurat_obj@misc)) {
    cli::cli_warn(
      "[{.fun GetVars2Regress}]: No QC columns found for regression in {.fun SCTransform}. Please run {.fun QCPatternDetect} first. Returning {.val NULL}."
    )

    return(NULL)
  }

  qc_cols <- unlist(seurat_obj@misc$qc_colnames, use.names = FALSE)
  meta <- seurat_obj[[]]

  existing_cols <- qc_cols[qc_cols %chin% colnames(meta)]

  if (length(existing_cols) == 0) {
    cli::cli_warn(
      "[{.fun GetVars2Regress}]: No matching QC columns found in seurat@meta.data Returning {.val NULL}."
    )

    return(NULL)
  }

  # 移除常量列（方差为0）
  valid_cols <- existing_cols[vapply(
    existing_cols,
    function(col) {
      x <- meta[[col]]
      !all(is.na(x)) && stats::var(x, na.rm = TRUE) > 0
    },
    logical(1L)
  )]

  if (length(valid_cols) == 0) {
    cli::cli_warn(
      "[{.fun GetVars2Regress}]: All QC columns are constant, skipping regression. Returning {.val NULL}."
    )
    return(NULL)
  }

  if (verbose && length(valid_cols) > 0) {
    cli::cli_text(
      "Using QC columns {valid_cols} for regression in {.fun SCTransform}"
    )
  }

  valid_cols
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
    cli::cli_alert_info(
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
        "[{.fun QCPatternDetect}]: Column {.val {.x}} already exists. Skipping pattern: {.val {.y}}"
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
    cli::cli_abort(c(
      "x" = "[{.fun QCFilter}]: seurat_obj must be a {.cls Seurat} object."
    ))
  }

  if (!is.list(data_filter.thresh)) {
    cli::cli_abort(c(
      "x" = "[{.fun QCFilter}]: {.arg data_filter.thresh} must be a named {.cls list}."
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
          cli::cli_warn(
            "[{.fun QCFilter}]: No threshold provided for {.val {col}}; skipping."
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
    cli::cli_abort(c(
      "x" = "[{.fun QCFilter}]: No valid filtering conditions generated."
    ))
  }

  full_expr <- purrr::reduce(
    all_conds,
    .f = function(x, y) rlang::expr(!!x & !!y)
  )

  # Evaluate safely
  logical_vec <- base::with(data = meta, expr = rlang::eval_tidy(full_expr))

  if (!is.logical(logical_vec) || length(logical_vec) != nrow(meta)) {
    cli::cli_abort(c(
      "x" = "[{.fun QCFilter}]: filtering condition did not produce a logical vector of length {.val {nrow(meta)}}."
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
#' @export
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
