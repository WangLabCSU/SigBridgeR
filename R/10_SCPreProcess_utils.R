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
