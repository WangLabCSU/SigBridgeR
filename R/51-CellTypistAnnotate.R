#' @title Annotate Cell Types Using CellTypist (Python Backend)
#' @description
#' Performs automated cell type annotation via the \code{celltypist} Python package
#' using \code{reticulate} integration. Accepts a \code{Seurat} object as input and
#' returns it enriched with CellTypist prediction results as metadata columns.
#'
#' Requires a Python environment with \code{celltypist} installed (see \url{https://github.com/Teichlab/celltypist}). The function
#' automatically attempts to locate a suitable Python interpreter, but users may
#' specify a custom path via the \code{python} argument.
#'
#' @param sc A \code{Seurat} object containing single-cell RNA-seq data. Must have
#'   a valid assay with gene expression matrix.
#' @param model Character. CellTypist model specification. One of:
#'   \itemize{
#'     \item Model name (e.g., \code{"Immune_All_Low"}): Loads the specified built-in model.
#'     \item Path to a local \code{.pkl} model file.
#'   }
#'
#' @param download Logical. Whether to automatically download the model first.
#'   Default: \code{TRUE}.
#' @param conda  Character. Conda environment name. Ignore `python` if specified.
#' @param python Character. Path to Python executable with \code{celltypist} installed.
#'   If \code{NULL} (default), auto-detected via \code{ListPyEnv()}. Must point to a valid Python binary.
#' @param venv_locations Character. Path to parent dirtectory storing Python virtual environment.
#' @param force_update Logical. Whether to force update the model file. `download` must be \code{TRUE} before this option is effective.
#' @param verbose Logical. Whether to print progress messages during annotation.
#'   Default: inherits from package option \code{getOption("SigBridgeR.verbose")}.
#' @param celltypist_tools Character. Path to the internal Python bridge script.
#'   Default: internal package resource (\code{system.file("python/73-CellTypistAnnotate.py", package = "SigBridgeR")}).
#'   Typically should not be modified by users.
#' @param ... Additional arguments passed to CellTypist's \code{annotate()} function via Python, such as:
#'   \itemize{
#'     \item \code{majority_voting}: Logical. Whether to refine predicted labels by running majority voting classifier after over-clustering. (Default: \code{FALSE})
#'     \item \code{mode}: Character. Prediction mode (\code{"best match"} or \code{"prob match"}). For \code{"best match"}, selects cell type with largest score; \code{"prob match"} enables multi-label classification. (Default: \code{"best match"})
#'     \item \code{p_thres}: Numeric. Probability threshold for multi-label classification in \code{"prob match"} mode. Ignored if \code{mode = "best match"}. (Default: 0.5)
#'     \item \code{transpose_input}: Logical. Whether to transpose input matrix. Set to \code{TRUE} if filename is in gene-by-cell format. (Default: \code{FALSE})
#'     \item \code{gene_file}: Character. Path to file with genes (one per line) corresponding to rows in mtx file. Ignored if input is not in mtx format.
#'     \item \code{cell_file}: Character. Path to file with cells (one per line) corresponding to columns in mtx file. Ignored if input is not in mtx format.
#'     \item \code{over_clustering}: Character or vector. Over-clustering specification: (1) path to plain file with one cluster ID per line; (2) metadata column name in AnnData; (3) vector/array of cluster assignments; or (4) omitted for heuristic approach. Ignored if \code{majority_voting = FALSE}.
#'     \item \code{use_GPU}: Logical. Whether to use GPU acceleration via rapids-singlecell for over-clustering. Only relevant when \code{majority_voting = TRUE}. (Default: \code{FALSE})
#'     \item \code{min_prop}: Numeric. Minimum proportion of dominant cell type required to name a subcluster. Subclusters below threshold are labeled 'Heterogeneous'. Ignored if \code{majority_voting = FALSE}. (Default: 0)
#'   }
#'
#' @return The input \code{Seurat} object with some metadata columns added, Column names may vary slightly depending on CellTypist version and options used. Usually cell type labels will be added to `meta.data`, and scoring matrix will be add to `misc$celltypist`
#'
#' @section Requirements:
#'   \itemize{
#'     \item R packages: \code{reticulate}, \code{AnnDataR}
#'     \item Python packages: \code{celltypist}, \code{scanpy}, \code{anndata}
#'     \item A working Python environment discoverable by \code{reticulate}
#'   }
#'
#'
#' @examples
#' \dontrun{
#' # Use a specific immune model with majority voting
#' annotated <- CellTypistAnnotate(
#'   seurat_obj,
#'   model = "Immune_All_Low",
#'   majority_voting = TRUE
#' )
#'
#' # Specify custom Python environment
#' annotated <- CellTypistAnnotate(
#'   seurat_obj,
#'   python = "/path/to/miniconda3/envs/celltypist/bin/python"
#' )
#' }
#' @family Single_Cell_Annotation_Method
#' @export
CellTypistAnnotate <- function(
  sc,
  model = NULL,
  download = TRUE,
  conda = NULL,
  python = NULL,
  venv_locations = NULL,
  force_update = TRUE,
  verbose = getFuncOption("verbose"),
  celltypist_tools = system.file(
    "python/73-CellTypistAnnotate.py",
    package = "SigBridgeR"
  ),
  ...
) {
  check_installed(c("anndataR", "reticulate"))
  chk::chk_is(sc, "Seurat")
  if (!is.null(model)) {
    chk::chk_character(model)
  }
  chk::chk_logical(download)
  chk::chk_file(celltypist_tools)

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green(
      "[CellTypist] Start annotating cell types"
    ))
  }

  py_envs <- ListPyEnv(venv_locations = venv_locations, verbose = FALSE)
  if (!is.null(conda)) {
    if (!is.null(python)) {
      cli::cli_warn("`python` is ignored due to `conda` specified")
    }
    SigBridgeRUtils::MatchArg(
      conda,
      py_envs$name[py_envs$type == "conda"],
      NULL
    )
    reticulate::use_condaenv(py_envs$python[py_envs$name == conda])
    if (verbose) {
      ts_cli$cli_alert_info("Using Conda Env: {.val {conda}}")
    }
  } else {
    python <- python %||% py_envs$python[[1]]
    chk::chk_file(python)
    reticulate::use_python(python)
    if (verbose) {
      ts_cli$cli_alert_info("Using Python: {.file {python}}")
    }
  }

  if (!is.null(model)) {
    if (!grepl("\\.pkl", model)) {
      model <- paste0(model, ".pkl")
    }
  }

  dots <- list2(...)

  # * initiate
  py <- reticulate::py
  reticulate::py_run_string("import os") # warm up

  py$adata <- anndataR::as_AnnData(
    x = sc,
    x_mapping = "data",
    output_class = "ReticulateAnnData"
  )
  py$model <- reticulate::r_to_py(model) # A string
  py$download <- reticulate::r_to_py(download) # A boolean
  py$verbose <- reticulate::r_to_py(verbose) # A boolean
  py$force_update <- reticulate::r_to_py(force_update) # A boolean
  if (length(dots) > 0) {
    var_names <- get_names_4_ids(..., .quoses = enquos(...))
    for (var_name in var_names) {
      py[[var_name]] <- reticulate::r_to_py(dots[[var_name]])
    }
  }
  reticulate::py_run_file(celltypist_tools)

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green("Annotation done"))
  }

  predicted_labels <- reticulate::py_to_r(py$predicted_labels)
  decision_matrix <- reticulate::py_to_r(py$decision_matrix)
  probability_matrix <- reticulate::py_to_r(py$probability_matrix)

  colnames(predicted_labels) <- paste0(
    "celltypist_",
    colnames(predicted_labels)
  )
  Seurat::AddMetaData(object = sc, metadata = predicted_labels) |>
    AddMisc(
      celltypist = list(
        decision_matrix = decision_matrix,
        probability_matrix = probability_matrix
      )
    )
}
