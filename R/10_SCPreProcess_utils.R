#' @title Compatible With Last Version
#' @family single_cell_preprocess
#' @keywords internal
compatible_with_3.0.2 <- function(..., params) {
  on.exit(gc(verbose = FALSE))
  dots <- rlang::list2(...)

  params$o$min.cells <- params$o$min.cells %||% dots$min_cells
  params$o$min.features <- params$o$min.cells %||% dots$min_features
  params$o$meta.data <- params$o$meta.data %||% dots$meta_data
  params$n$normalization.method <- params$n$normalization.method %||%
    dots$normalization_method
  params$n$scale.factor <- params$n$scale.factor %||% dots$scale_factor
  params$s$scale.features <- params$s$scale.features %||% dots$scale_features
  params$v$selection.method <- params$v$selection.method %||%
    dots$selection_method
  params$c$resolution <- params$c$resolution %||% dots$resolution

  params
}

#' @keywords internal
#' @family single_cell_preprocess
sc_processing_strategies <- list(
  o = function(object, params) {
    rlang::exec(
      SeuratObject::CreateSeuratObject,
      counts = object,
      !!!params
    )
  },
  n = function(object, params) {
    rlang::exec(
      Seurat::NormalizeData,
      object = object,
      !!!params
    )
  },
  v = function(object, params) {
    rlang::exec(
      Seurat::FindVariableFeatures,
      object = object,

      !!!params
    )
  },
  s = function(object, params) {
    rlang::exec(
      Seurat::ScaleData,
      object = object,

      !!!params
    )
  },
  p = function(object, params) {
    rlang::exec(
      Seurat::RunPCA,
      object = object,

      !!!params
    )
  },
  c = function(object, params) {
    rlang::exec(
      Seurat::FindClusters,
      object = object,

      !!!params
    )
  },
  e = function(object, params) {
    rlang::exec(
      Seurat::FindNeighbors,
      object = object,
      !!!params
    )
  },
  t = function(object, params) {
    rlang::exec(
      Seurat::RunTSNE,
      object = object,
      !!!params
    )
  },
  u = function(object, params) {
    rlang::exec(
      Seurat::RunUMAP,
      object = object,
      !!!params
    )
  },
  r = function(object, params) {
    rlang::exec(
      Seurat::SCTransform,
      object = object,
      !!!params
    )
  },
  i = function(object, params) {
    rlang::exec(
      Seurat::IntegrateLayers,
      object = object,
      !!!params
    )
  }
)


#' @keywords internal
#' @family single_cell_preprocess
has_pattern <- function(qc_list) {
  !is.null(qc_list) && !is.null(qc_list$pattern) && length(qc_list$pattern) > 0
}

#' @keywords internal
#' @family single_cell_preprocess
is_filtering <- function(filter_list) {
  !is.null(filter_list)
}

#' @keywords internal
#' @family single_cell_preprocess
modifyList2 <- function(x, val) {
  chk::chk_list(x)
  chk::chk_list(val)

  purrr::imap(val, function(v_val, key) {
    if (key %chin% names(x) && is.list(x[[key]]) && is.list(v_val)) {
      return(utils::modifyList(x[[key]], v_val))
    }
    v_val
  })
}

#' @keywords internal
#' @family single_cell_preprocess
handle_usr_params <- function(usr_params) {
  default <- list(
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
  )
  modifyList2(default, usr_params)
}
