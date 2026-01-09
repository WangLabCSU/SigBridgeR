#' @title Integrate Single-Cell Datasets
#'
#' @description
#' Unified interface for integrating heterogeneous single-cell datasets.
#' Supports \code{Seurat} objects and matrices (via gene-wise aggregation and column concatenation).
#'
#' @param ... For matrix method: matrices to integrate, optionally named.
#'   For Seurat method: Seurat objects to integrate, followed by additional
#'   parameters passed to integration functions, Seurat objects and parameters are automatically filtered and passed to appropriate functions.
#' @param pipeline Character string specifying processing steps for Seurat integration.
#'   Each letter represents a step: \code{"n"} = NormalizeData, \code{"s"} = ScaleData,
#'   \code{"v"} = FindVariableFeatures, \code{"p"} = RunPCA, \code{"i"} = IntegrateLayers,
#'   \code{"e"} = FindNeighbors (because `n` was used),
#'   \code{"c"} = FindClusters, \code{"t"} = RunTSNE, \code{"u"} = RunUMAP,
#'   \code{"r"} = SCTransform (because `t` was used),. Default: \code{"nsvpiectu"}.
#' @param add.cell.ids Character vector of prefixes for cell IDs when merging
#'   Seurat objects. Auto-inferred from argument names if not provided. (See Examples for details)
#' @param collapse If TRUE, merge layers of the same name together; if FALSE, appends labels to the layer name
#' @param merge.data Merge the data slots instead of just merging the counts (which requires renormalization); this is recommended if the same normalization approach was applied to all objects
#' @param merge.dr Choose how to handle merging dimensional reductions: - “TRUE”: merge dimensional reductions with the same name across objects; dimensional reductions with different names are added as-is - “NA”: keep dimensional reductions from separate objects separate; will append the project name for duplicate reduction names - “FALSE”: do not add dimensional reductions
#' @param project Project name for the Seurat object
#' @param .quos Ignore it. Please do not pass any value
#'
#'
#' @return
#' \itemize{
#'   \item \strong{Matrix method}: A matrix with all genes (union) and
#'     concatenated cells. Missing values filled with \code{NA}.
#'   \item \strong{Seurat method}: Integrated \code{Seurat} object with
#'     processed layers according to \code{pipeline}.
#' }
#'
#'
#' @examples
#' \dontrun{
#' # Matrix integration
#' mat1 <- matrix(rpois(100, 5), nrow = 20,
#'                dimnames = list(paste0("G", 1:20), paste0("C", 1:5)))
#' mat2 <- matrix(rpois(120, 6), nrow = 20,
#'                dimnames = list(paste0("G", 11:30), paste0("C", 1:6)))
#' integrated <- SCIntegrate(mat1, mat2)
#' dim(integrated)  # 20 genes × 11 cells
#'
#' # Add customed prefixes to colnames
#' integrated <- SCIntegrate(A = mat1, B = mat2)
#' # > colnames(integrated)
#' # [1] "A_C1" "A_C2" "A_C3" "A_C4" "A_C5" "B_C1" "B_C2" "B_C3" "B_C4" "B_C5" "B_C6"
#'
#' # Seurat integration with custom pipeline
#' seu1 <- CreateSeuratObject(mat1)
#' seu2 <- CreateSeuratObject(mat2)
#' integrated_seu <- SCIntegrate(seu1, seu2, pipeline = "nsfpi")
#' }
#'
#' @name SCIntegrate
#' @export
NULL

#' @rdname SCIntegrate
#' @export
SCIntegrate <- function(
  ...,
  pipeline = "nsvpiectu",
  add.cell.ids = NULL,
  collapse = FALSE,
  merge.data = TRUE,
  merge.dr = FALSE,
  project = "Integrated Seurat"
) {
  # ! don't use S3 method
  dots <- rlang::list2(...)
  if (length(dots) == 0) {
    cli::cli_abort(c("[{.fun SCIntegrate}]: No arguments provided."))
  }
  .quos <- rlang::enquos(...)
  if (inherits(dots[[1]], "matrix")) {
    return(SCIntegrate.matrix(..., .quos = .quos))
  }
  if (inherits(dots[[1]], "Seurat")) {
    return(SCIntegrate.Seurat(
      ...,
      pipeline = pipeline,
      add.cell.ids = add.cell.ids,
      collapse = collapse,
      merge.data = merge.data,
      merge.dr = merge.dr,
      project = project,
      .quos = .quos
    ))
  }

  cls <- c("Seurat", "matrix")
  cli::cli_abort(c(
    "x" = "[{.fun SCIntegrate}]: No implementation for class {.cls {class(dots[[1]])}}",
    ">" = "Available classes: {.cls {cls}}"
  ))
}

#' @rdname SCIntegrate
#' @section Matrix Method Details:
#'   \itemize{
#'     \item Duplicate genes (e.g., due to symbol aliasing) are resolved via
#'           \code{\link{AggregateDups}} (default: sum).
#'     \item Dataset identifiers are auto-inferred from argument names or
#'           object names (e.g., \code{SCIntegrate(matA, B = matB)} → prefixes \code{"matA"}, \code{"B"}).
#'   }
#' @examples
#' \dontrun{
#' mat1 <- matrix(rpois(100, 5), nrow = 20, dimnames = list(paste0("G", 1:20), paste0("C", 1:5)))
#' mat2 <- matrix(rpois(120, 6), nrow = 20, dimnames = list(paste0("G", 1:20), paste0("C", 1:6)))
#' integrated_mat <- SCIntegrate(mat1, mat2)
#' dim(integrated_mat)  # 20 genes × 11 cells
#' head(colnames(integrated_mat))
#' }
#' @export
SCIntegrate.matrix <- function(..., .quos = NULL) {
  dots <- rlang::list2(...)
  .quos <- .quos %||% rlang::enquos(...)
  is_mat <- vapply(dots, \(x) inherits(x, "matrix"), logical(1))
  # * remove dups
  mats <- purrr::map(dots[is_mat], \(x) {
    rlang::exec(
      AggregateDups,
      x = x,
      verbose = FALSE,
      !!!SigBridgeRUtils::FilterArgs4Func(dots[!is_mat], AggregateDups)
    )
  })

  all_genes <- unique(unlist(lapply(mats, rownames)))

  prefixes <- get_names_4_ids(..., .quoses = .quos)[is_mat]

  # * merge
  dt_list <- vector("list", length(mats))
  for (i in seq_along(mats)) {
    mat <- mats[[i]]
    prefix <- prefixes[i]

    # 转换为data.table
    dt <- data.table::as.data.table(mat, keep.rownames = "gene")

    # 给列名添加前缀(gene)
    col_names <- names(dt)[-1]
    data.table::setnames(dt, col_names, paste0(prefix, "_", col_names))

    dt_list[[i]] <- dt
  }

  result_dt <- Reduce(
    function(x, y) {
      merge(x, y, by = "gene", all = TRUE)
    },
    dt_list
  )
  result_mat <- as.matrix(result_dt[, -1])
  rownames(result_mat) <- result_dt$gene

  result_mat
}

#' @rdname SCIntegrate
#' @export
SCIntegrate.Seurat <- function(
  ...,
  pipeline = "nsvpiectu",
  add.cell.ids = NULL,
  collapse = FALSE,
  merge.data = TRUE,
  merge.dr = FALSE,
  project = "Integrated Seurat",
  .quos = NULL
) {
  chk::chk_character(pipeline)
  chk::chk_length(pipeline)
  # * parameters
  dots <- rlang::list2(...)
  .quos <- .quos %||% rlang::enquos(...)
  method <- dots$method %||% Seurat::HarmonyIntegration
  verbose <- dots$verbose %||%
    SigBridgeRUtils::getFuncOption("verbose") %||%
    TRUE

  is_seurat <- vapply(dots, \(x) inherits(x, "Seurat"), logical(1))

  # * merge
  if (verbose) {
    cli::cli_text("Start merging Seurat objects")
  }
  first_seurat <- dots[[1]]
  other_seurat <- unlist(dots[is_seurat][-1])
  merged <- merge(
    x = first_seurat,
    y = other_seurat,
    add.cell.ids = add.cell.ids %||%
      get_names_4_ids(..., .quoses = .quos)[is_seurat],
    collapse = collapse,
    merge.data = merge.data,
    merge.dr = merge.dr,
    project = project
  )

  if (!is.function(method)) {
    # ? handling user-provided method
    cli::cli_warn(
      "[{.fun SCIntegrate.Seurat}]: {.arg method} must be a function, returning merged Seurat object without integrated layers"
    )
    return(merged)
  }

  # * use `pipeline` to control steps
  step_map <- rlang::list2(
    r = function() {
      rlang::exec(
        Seurat::SCTransform,
        object = merged,
        !!!SigBridgeRUtils::FilterArgs4Func(
          dots[!is_seurat],
          Seurat::SCTransform
        )
      )
    },
    n = function() {
      rlang::exec(
        Seurat::NormalizeData,
        object = merged,
        !!!SigBridgeRUtils::FilterArgs4Func(
          dots[!is_seurat],
          Seurat::NormalizeData
        )
      )
    },
    s = function() {
      rlang::exec(
        Seurat::ScaleData,
        object = merged,
        !!!SigBridgeRUtils::FilterArgs4Func(
          dots[!is_seurat],
          Seurat::ScaleData
        )
      )
    },
    v = function() {
      rlang::exec(
        Seurat::FindVariableFeatures,
        object = merged,
        !!!SigBridgeRUtils::FilterArgs4Func(
          dots[!is_seurat],
          method
        )
      )
    },
    p = function() {
      rlang::exec(
        Seurat::RunPCA,
        object = merged,
        !!!SigBridgeRUtils::FilterArgs4Func(
          dots[!is_seurat],
          Seurat::RunPCA
        )
      )
    },
    i = function() {
      rlang::exec(
        Seurat::IntegrateLayers,
        object = merged,
        !!!dots[!is_seurat] # pass ... param to method
      )
    },
    e = function() {
      rlang::exec(
        FindNeighbors,
        object = merged,

        !!!SigBridgeRUtils::FilterArgs4Func(
          dots[!is_seurat],
          FindNeighbors
        )
      )
    },
    c = function() {
      rlang::exec(
        Seurat::FindClusters,
        object = merged,
        !!!SigBridgeRUtils::FilterArgs4Func(
          dots[!is_seurat],
          Seurat::FindClusters
        )
      )
    },
    t = function() {
      rlang::exec(
        Seurat::RunTSNE,
        object = merged,
        !!!SigBridgeRUtils::FilterArgs4Func(
          dots[!is_seurat],
          Seurat::RunTSNE
        )
      )
    },
    u = function() {
      rlang::exec(
        Seurat::RunUMAP,
        object = merged,
        !!!SigBridgeRUtils::FilterArgs4Func(
          dots[!is_seurat],
          Seurat::RunUMAP
        )
      )
    }
  )

  steps_to_run <- unlist(strsplit(pipeline, ""))
  execution_queue <- step_map[steps_to_run[steps_to_run %chin% names(step_map)]]

  for (step_name in names(execution_queue)) {
    merged <- execution_queue[[step_name]]()
  }

  if (verbose) {
    cli::cli_alert_success("Integration finished")
  }

  merged
}

#' @title Get Names for Dot Variables
#' @description Extract or infer names from \code{...} arguments.
#' @param ... Arguments to extract names from.
#' @param .quoses predefined
#' @return Character vector of names.
#' @keywords internal
#' @examples
#' \dontrun{
#' get_names_4_ids(a = 1, b = 2)
#' # [1] "a" "b"
#' c <- 3
#' get_names_4_ids(a = 1, b = 2, c)
#' # [1] "a" "b" "c"
#' get_names_4_ids()
#' # character(0)
#'
#' mat1 <- matrix(1:2)
#' mat2 <- matrix(3:4)
#' get_names_4_ids(mat1, mat2)
#' # [1] "mat1" "mat2"
#' get_names_4_ids(mat1, c = mat2)
#' # [1] "mat1" "c"
#' }
get_names_4_ids <- function(..., .quoses = NULL) {
  dots <- rlang::list2(...)

  # 如果提供了预计算的 quosures 就用它，否则当场捕获
  if (is.null(.quoses)) {
    quoses <- rlang::enquos(..., .named = TRUE)
  } else {
    quoses <- .quoses
  }

  var_names <- names(dots)
  if (is.null(var_names)) {
    var_names <- rep("", length(dots))
  }
  unnamed <- which(var_names == "")

  if (length(unnamed) > 0) {
    var_names[unnamed] <- purrr::map_chr(quoses[unnamed], rlang::as_label)
  }

  var_names
}
