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
has_pattern <- function(qc_list) {
  if (!is.list(qc_list)) {
    return(FALSE)
  }
  !is.null(qc_list) && !is.null(qc_list$pattern) && length(qc_list$pattern) > 0
}

#' @keywords internal
#' @family single_cell_preprocess
is_filtering <- function(filter_list) {
  if (!is.list(filter_list)) {
    return(FALSE)
  }
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
