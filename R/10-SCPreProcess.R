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
#' Default: \code{assay} ("RNA"), \code{nFeature_RNA} (200-6000), \code{nCount_RNA} (500-50000), \code{percent.mt} (<20), \code{percent.rp} (<60).
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
#' @family input_preprocess
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
  sc = NULL,
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
    # assay = "RNA", # detect default assay
    nFeature_thresh = c(200L, 6000L),
    nCount_thresh = c(500L, 50000L),
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
  dots <- list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose") %||% TRUE

  params <- handle_usr_params(params) # combine with default params

  # * handling pipeline
  steps <- unlist(strsplit(pipeline, ""))
  steps_to_run <- steps[steps %chin% names(SCPreProcessStrategy)] # a vector

  unknown <- setdiff(steps, steps_to_run)
  if (length(unknown) != 0) {
    cli::cli_abort(c(
      "x" = "[{.fun SCPreProcess}]: Unknown pipeline steps: {.val {unknown}}"
    ))
  }
  if (steps_to_run[[1]] == "o" && !is.null(sc)) {
    sc_seurat <- if (!is.null(params$o$counts)) {
      cli::cli_warn(
        "[{.fun SCPreProcess}]: The parameter `params$o$counts` is not NULL, `sc` will be ignored."
      )
      exec(
        SeuratObject::CreateSeuratObject,
        !!!params$o
      )
    } else {
      exec(
        SeuratObject::CreateSeuratObject,
        counts = sc,
        !!!params$o
      )
    }
  } else {
    # * CreateSeuratObject("o") is not in the pipeline. Use the provided custom method
    if (verbose) {
      cli::cli_inform("Use custom method to create Seurat object")
    }
    letter <- steps_to_run[[1]]

    if (!is.null(sc)) {
      cli::cli_warn(
        "It is recommend to set `sc` to NULL when using custom method"
      )
    }

    sc_seurat <- exec(
      SCPreProcessStrategy[[letter]],
      sc,
      !!!params[[letter]]
    )
  }

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
      data_filter$nFeature_thresh[1],
      data_filter$nFeature_thresh[2]
    )
    chk::chk_lt(
      data_filter$nCount_thresh[1],
      data_filter$nCount_thresh[2]
    )
    sc_seurat <- QCFilter(
      seurat_obj = sc_seurat,
      data_filter.thresh = data_filter,
      verbose = verbose
    )
  }

  # the first step is CreateSeuratObject("o")
  for (step in 2:length(steps_to_run)) {
    # step: an index
    letter <- steps_to_run[[step]]
    step_fun <- SCPreProcessStrategy[[letter]] # function

    # if (!is.function(step_fun)) {
    #   cli::cli_abort(c(
    #     "x" = "[{.fun SCPreProcess}]: Step {.val {steps_to_run[[step]]}} function not found",
    #     ">" = "Available steps: {.val {names(SCPreProcessStrategy)}}",
    #     ">" = "Use {.fun RegisterSeuratMethod} to add custom function"
    #   ))
    # }
    sc_seurat <- exec(step_fun, sc_seurat, !!!params[[letter]])
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
    # assay = "RNA", # detect default assay
    nFeature_thresh = c(200L, 6000L),
    nCount_thresh = c(500L, 50000L),
    percent.mt = 20L, # mitochondrial genes
    percent.rp = 60L # ribosomal protein genes
  ),
  column2only_tumor = NULL
) {
  check_installed("dplyr")
  # Both `anndata` and `anndataR` are based on R6

  # * check
  if (!is.null(column2only_tumor)) {
    chk::chk_character(column2only_tumor)
  }
  chk::chk_character(pipeline)
  purrr::walk(c(quality_control, data_filter, params), ~ chk::chk_list)
  purrr::walk(params, ~ chk::chk_list)

  # * dots arguments, to be compatible last version
  dots <- list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose") %||% TRUE

  params <- handle_usr_params(params) # combine with default params

  # * handling pipeline
  steps <- unlist(strsplit(pipeline, ""))
  steps_to_run <- steps[steps %chin% names(SCPreProcessStrategy)] # a vector

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
    if (!is.null(params$o$meta.data)) {
      chk::chk_data(params$o$meta.data)
      seurat <- SeuratObject::AddMetaData(seurat, params$o$meta.data)
    }
    if (params$o$min.cells != 0) {
      chk::chk_whole_number(params$o$min.cells)
      gene_cell_counts <- SigBridgeRUtils::rowSums3(
        SeuratObject::LayerData(seurat, layer = "counts") > 0
      )
      seurat <- seurat[gene_cell_counts >= params$o$min.cells, ]
    }
    if (params$o$min.features != 0) {
      chk::chk_whole_number(params$o$min.features)
      cell_gene_counts <- SigBridgeRUtils::colSums3(
        SeuratObject::LayerData(seurat, layer = "counts") > 0
      )
      seurat <- seurat[, cell_gene_counts >= params$o$min.features]
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

    seurat <- exec(
      SeuratObject::CreateSeuratObject,
      counts = count_mat,
      !!!params$o
    )

    if (ncol(sc$var) != 0) {
      seurat <- AddMetaFeature(seurat, sc$var)
    }
    seurat
  } else {
    Abort(
      "x" = "{.arg sc} must be an anndata or anndataR object",
      ">" = "Current input is of class {.cls {class(sc)}}"
    )
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
      data_filter$nFeature_thresh[1],
      data_filter$nFeature_thresh[2]
    )
    chk::chk_lt(
      data_filter$nCount_thresh[1],
      data_filter$nCount_thresh[2]
    )
    sc_seurat <- QCFilter(
      seurat_obj = sc_seurat,
      data_filter.thresh = data_filter,
      verbose = verbose
    )
  }

  # a list containing functions
  execution_queue <- SCPreProcessStrategy[steps_to_run]

  # the first stop is CreateSeuratObject("o")
  # step: an index
  for (step in 2:length(steps_to_run)) {
    # step_fun: a function
    step_fun <- execution_queue[[step]]
    sc_seurat <- exec(
      step_fun,
      sc_seurat,
      # steps_to_run[[step]]: a letter
      !!!params[[steps_to_run[[step]]]]
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
  dots <- list2(...)
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
