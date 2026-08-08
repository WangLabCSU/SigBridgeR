#' @title Annotate Cell Types Using CellTypist (Python Backend)
#' @family Single_Cell_Annotation_Method
#'
#' @description
#' Performs automated cell type annotation via the \code{celltypist} Python package
#' using \code{reticulate} integration. Accepts a \code{Seurat} object as input and
#' returns it enriched with CellTypist prediction results as metadata columns.
#'
#' Requires a Python environment with \code{celltypist} installed (see
#' \url{https://github.com/Teichlab/celltypist}). By default the function locates a
#' suitable Python interpreter automatically via \code{reticulate}; use
#' \code{reticulate::use_python()} to pin a specific environment beforehand.
#'
#' @param sc A \code{Seurat} object containing single-cell RNA-seq data. Must have
#'   a valid assay with gene expression matrix.
#' @param model Character. CellTypist model specification. One of:
#'   * A built-in model name (e.g., `"Immune_All_Low"`); a `.pkl` suffix is appended
#'     automatically if missing.
#'   * A path to a local `.pkl` model file.
#' @param download Logical. Whether to automatically download the model first.
#'   Default: `TRUE`.
#' @param force_update Logical. Whether to force re-downloading the model files
#'   when \code{download = TRUE}. Default: `FALSE`.
#' @param majority_voting Logical. Whether to refine predicted labels by running a
#'   majority voting classifier after over-clustering. Default: `FALSE`.
#' @param mode Character. Prediction mode: `"best match"` selects the cell type with
#'   the largest score; `"prob match"` enables multi-label classification.
#'   Default: `"best match"`.
#' @param p_thres Numeric. Probability threshold for multi-label classification in
#'   `"prob match"` mode. Ignored when \code{mode = "best match"}. Default: `0.5`.
#' @param transpose_input Logical. Whether to transpose the input matrix. Set to
#'   `TRUE` if the input is in gene-by-cell format. Default: `FALSE`.
#' @param gene_file Character. Path to a file with genes (one per line) corresponding
#'   to rows of the input matrix. Default: `NULL`.
#' @param cell_file Character. Path to a file with cells (one per line) corresponding
#'   to columns of the input matrix. Default: `NULL`.
#' @param over_clustering Character or vector. Over-clustering specification used for
#'   majority voting: (1) path to a plain file with one cluster ID per line;
#'   (2) metadata column name in AnnData; (3) vector/array of cluster assignments;
#'   or (4) omitted for the heuristic approach. Ignored when
#'   \code{majority_voting = FALSE}. Default: `NULL`.
#' @param use_GPU Logical. Whether to use GPU acceleration (via rapids-singlecell)
#'   for over-clustering. Only relevant when \code{majority_voting = TRUE}.
#'   Default: `FALSE`.
#' @param min_prop Numeric. Minimum proportion of the dominant cell type required to
#'   name a subcluster; subclusters below the threshold are labeled `'Heterogeneous'`.
#'   Ignored when \code{majority_voting = FALSE}. Default: `0L`.
#' @param verbose Logical. Whether to print progress messages during annotation.
#'   Default: \code{getFuncOption("verbose")}.
#' @param assay which Seurat assay to use
#' @param ... Additional arguments passed through to the CellTypist Python backend.
#'
#' @return The input \code{Seurat} object with prediction results added:
#'   cell type labels are stored as new columns in `meta.data` (prefixed with
#'   `celltypist_`), and the scoring matrices are stored in `misc$celltypist`.
#'   Column names may vary slightly depending on the CellTypist version and options used.
#'
#' @section Requirements:
#'
#' * R packages: `reticulate`, `AnnDataR`
#' * Python packages: `celltypist`, `scanpy`, `anndata`
#' * A working Python environment discoverable by `reticulate`
#'
#' @examples
#' \dontrun{
#' # Use a built-in immune model with majority voting
#' annotated <- CellTypistAnnotate(
#'   sc = seurat_obj,
#'   model = "Immune_All_Low",
#'   majority_voting = TRUE
#' )
#'
#' # Probability-based multi-label annotation
#' annotated <- CellTypistAnnotate(
#'   sc = seurat_obj,
#'   model = "Immune_All_Low",
#'   mode = "prob match",
#'   p_thres = 0.5
#' )
#'
#' # Use a local model file
#' annotated <- CellTypistAnnotate(
#'   sc = seurat_obj,
#'   model = "/path/to/My_Model.pkl"
#' )
#' }
#' @export
CellTypistAnnotate <- function(
  sc,
  model = NULL,
  download = TRUE,
  force_update = FALSE,
  majority_voting = FALSE,
  mode = c("best match", "prob match"),
  p_thres = 0.5,
  transpose_input = FALSE,
  gene_file = NULL,
  cell_file = NULL,
  over_clustering = NULL,
  use_GPU = FALSE,
  min_prop = 0L,
  verbose = getFuncOption("verbose"),
  assay = "RNA",
  ...
) {
  check_installed("anndataR")
  chk::chk_is(sc, "Seurat")
  if (!is.null(model)) {
    chk::chk_character(model)
    if (!endsWith(model, ".pkl")) {
      model <- paste0(model, ".pkl")
    }
  }
  mode <- arg_match(mode)
  chk::chk_logical(download)
  chk::chk_logical(force_update)
  chk::chk_logical(majority_voting)
  chk::chk_numeric(p_thres)
  chk::chk_logical(transpose_input)
  chk::chk_logical(use_GPU)
  chk::chk_numeric(min_prop)
  chk::chk_logical(verbose)
  if (!is.null(gene_file)) {
    chk::chk_character(gene_file)
  }
  if (!is.null(cell_file)) {
    chk::chk_character(cell_file)
  }
  # ----------------------------------------------------------------------------

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green(
      "[CellTypist] Start annotating cell types"
    ))
  }

  anno_res <- PyModule$celltypist$annotate_celltypist(
    adata = anndataR::as_AnnData(
      x = sc,
      x_mapping = "data",
      assay_name = assay,
      output_class = "ReticulateAnnData"
    ),
    model = model,
    ...,
    transpose_input = transpose_input,
    gene_file = gene_file,
    cell_file = cell_file,
    majority_voting = majority_voting,
    over_clustering = over_clustering,
    use_GPU = use_GPU,
    mode = mode,
    p_thres = p_thres,
    min_prop = min_prop,
    download = download,
    force_update = force_update,
    verbose = verbose
  ) |>
    reticulate::py_to_r()

  predicted_labels <- anno_res[[1L]]
  decision_matrix <- anno_res[[2L]]
  probability_matrix <- anno_res[[3L]]

  colnames(predicted_labels) <- paste0(
    "celltypist_",
    colnames(predicted_labels)
  )
  Seurat::AddMetaData(object = sc, metadata = predicted_labels) |>
    AddMisc(
      celltypist = list(
        decision_matrix = decision_matrix,
        probability_matrix = probability_matrix
      ),
      cover = FALSE
    )
}

#' @title Download CellTypist Models
#' @family Single_Cell_Annotation_Method
#'
#' @description
#' Downloads CellTypist pre-trained models via the \pkg{celltypist} Python
#' backend. A single model, multiple models, or all models can be downloaded.
#' The models are saved to CellTypist's local models directory and can
#' subsequently be used by \code{\link{CellTypistAnnotate}} without
#' re-downloading each time.
#'
#' To see all available models and their descriptions, see
#' \url{https://github.com/Teichlab/celltypist}.
#'
#' @param model Character. Name(s) of the CellTypist model(s) to download,
#'   e.g. `"Immune_All_Low"`. A `.pkl` suffix is appended automatically if
#'   missing. When `NULL`, all models are downloaded. Default: `NULL`.
#' @param force_update Logical. Whether to force re-downloading the CellTypist
#'   model files even when they already exist locally. Default: `FALSE`.
#'
#' @return Invisibly `NULL`. The models are downloaded to the CellTypist models
#'   directory (see `models.models_path` in the Python backend).
#'
#' @section Requirements:
#'
#' * A Python environment with `celltypist` installed, discoverable by `reticulate`.
#'
#' @examples
#' \dontrun{
#' # Download all CellTypist models
#' DownloadCellTypistModels()
#'
#' # Download a single model
#' DownloadCellTypistModels(model = "Immune_All_Low")
#'
#' # Force re-downloading a specific model
#' DownloadCellTypistModels(model = "Immune_All_Low", force_update = TRUE)
#' }
#' @export
DownloadCellTypistModels <- function(model = NULL, force_update = FALSE) {
  chk::chk_logical(force_update)
  if (!is.null(model)) {
    chk::chk_character(model)
    model <- ifelse(
      endsWith(model, ".pkl"),
      model,
      paste0(model, ".pkl")
    )
  }
  PyModule$celltypist$download_models(
    model = model,
    force_update = force_update
  )
}

#' @title List Available CellTypist Models
#' @family Single_Cell_Annotation_Method
#'
#' @description
#' Lists the CellTypist pre-trained models available locally. By default it
#' reads all model files under `~/.celltypist/data/models`, the default
#' CellTypist models directory (overridable via the `CELLTYPIST_FOLDER`
#' environment variable). Models are downloaded with
#' \code{\link{DownloadCellTypistModels}} and used by
#' \code{\link{CellTypistAnnotate}}.
#'
#' @param path Character. Path to the CellTypist models directory. Default:
#'   `NULL`, which resolves to `$CELLTYPIST_FOLDER/data/models` (or
#'   `~/.celltypist/data/models` when the environment variable is unset).
#' @param verbose Logical. Whether to print the model list to the console.
#'   Default: \code{getFuncOption("verbose")}.
#'
#' @return A character vector of available model names (without the `.pkl`
#'   suffix). Invisibly returned when \code{verbose = TRUE}.
#'
#' @examples
#' \dontrun{
#' # List all locally available CellTypist models
#' ListCellTypistModels()
#'
#' # List models in a custom directory, without printing
#' ListCellTypistModels(path = "~/my_models", verbose = FALSE)
#' }
#' @export
ListCellTypistModels <- function(
  path = NULL,
  verbose = getFuncOption("verbose")
) {
  chk::chk_logical(verbose)
  path <- path %||%
    file.path(
      Sys.getenv(
        "CELLTYPIST_FOLDER",
        unset = path.expand("~/.celltypist")
      ),
      "data",
      "models"
    )
  chk::chk_string(path)

  if (!dir.exists(path)) {
    Abort(
      "CellTypist models directory not found: {.path {path}}.",
      "Run {.fn DownloadCellTypistModels} to download the models first."
    )
  }

  models <- sort(list.files(
    path = path,
    pattern = "[.]pkl$"
  ))
  model_names <- gsub("[.]pkl$", "", models)

  if (verbose) {
    if (length(model_names)) {
      cli::cli_h2("Available CellTypist models ({length(model_names)})")
      cli::cli_li(model_names)
    } else {
      cli::cli_alert_info(
        "No CellTypist models found in {.path {path}}."
      )
    }
    invisible(model_names)
  } else {
    model_names
  }
}
