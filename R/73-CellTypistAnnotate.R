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
#'     \item \code{NULL} (default): Uses CellTypist's default model (auto-downloaded if \code{download = TRUE}).
#'     \item Model name (e.g., \code{"Immune_All_Low"}): Downloads or loads the specified built-in model.
#'     \item Path to a local \code{.pkl} model file.
#'   }
#'   See \url{https://celltypist.readthedocs.io/en/latest/models.html} for available models.
#' @param download Logical. Whether to automatically download the specified model first.
#'   Default: \code{TRUE}.
#' @param conda  Character. Conda environment name. Ignore `python` if specified.
#' @param python Character. Path to Python executable with \code{celltypist} installed.
#'   If \code{NULL} (default), auto-detected via \code{ListPyEnv()}. Must point to a valid Python binary.
#' @param verbose Logical. Whether to print progress messages during annotation.
#'   Default: inherits from package option \code{getOption("SigBridgeR.verbose")}.
#' @param celltypist_tools Character. Path to the internal Python bridge script.
#'   Default: internal package resource (\code{system.file("python/73-CellTypistAnnotate.py", package = "SigBridgeR")}).
#'   Typically should not be modified by users.
#' @param ... Additional arguments passed to CellTypist's \code{annotate()} function via Python, such as:
#'   \itemize{
#'     \item \code{majority_voting}: Logical. Whether to apply majority voting for cluster-level annotation.
#'     \item \code{mode}: Character. Prediction mode (\code{"best match"} or \code{"prob match"}).
#'     \item \code{p_thres}: Numeric. Probability threshold for \code{"prob match"} mode.
#'   }
#'   See \url{https://celltypist.readthedocs.io/en/latest/api.html#celltypist.annotate} for full options.
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
  verbose = getFuncOption("verbose"),
  celltypist_tools = system.file(
    "python/73-CellTypistAnnotate.py",
    package = "SigBridgeR"
  ),
  ...
) {
  rlang::check_installed(c("anndataR", "reticulate"))
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

  py_envs <- ListPyEnv(verbose = FALSE)
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

  dots <- rlang::list2(...)

  # * initiate
  py <- reticulate::py

  py$adata <- anndataR::as_AnnData(
    x = sc,
    x_mapping = "data",
    output_class = "ReticulateAnnData"
  )
  py$model <- reticulate::r_to_py(model) # A string
  py$download <- reticulate::r_to_py(download)
  py$verbose <- reticulate::r_to_py(verbose) # A boolean
  if (length(dots) > 0) {
    var_names <- get_names_4_ids(..., .quoses = rlang::enquos(...))
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
  Seurat::AddMetaData(object = sc, metadata = predicted_labels) %>%
    AddMisc(
      celltypist = list(
        decision_matrix = decision_matrix,
        probability_matrix = probability_matrix
      )
    )
}
