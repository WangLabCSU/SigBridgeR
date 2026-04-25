#' @title Screen Function Template
#'
#' @param matched_bulk Matrix or data frame of preprocessed bulk RNA-seq expression
#'        data (genes x samples). Column names must match names/IDs in `phenotype`.
#' @param sc_data A matrix/Matrix (genes x cells) or a Seurat object containing scRNA-seq data to be screened.
#' @param phenotype Phenotype data, either:
#'        - Named vector (names match `matched_bulk` columns), or
#'        - Patient survival Data frame with row names matching `matched_bulk` columns, colnames named "time" and "status"
#' @param label_type Character specifying phenotype label type
#' @param phenotype_class Type of phenotypic outcome (must be consistent with input data):
#'        - `"binary"`: Binary traits (e.g., case/control)
#'        - `"continuous"`: Continuous measurements
#'        - `"survival"`: Survival infomation
#' @param ... Additional arguments passed to the function. Common parameters include:
#'   \describe{
#'     \item{verbose}{Logical. Whether to print verbose output (default: TRUE).}
#'   }
#'
#' @return A named list containing:
#'   \describe{
#'     \item{scRNA_data}{Modified single-cell data object with integrated screening results.}
#'   }
#' @family scSurvival
#' @export
DoscSurvival <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = NULL,
  phenotype_class = c("survival"),
  env_params = list(),
  scsurvival_params = list(),
  ...
) {
  CheckInstalled("Exceret/rscSurvival")
  chk::chk_is(sc_data, "Seurat")
  chk::chk_character(label_type)
  phenotype_class <- SigBridgeRUtils::MatchArg(
    phenotype_class,
    c("survival"),
    NULL
  )
  chk::chk_number(maxiter)
  chk::chk_number(tred)
  # scAB can't tolerate NA
  chk::chk_not_any_na(matched_bulk)
  chk::chk_not_any_na(phenotype)

  # Extract additional arguments
  dots <- rlang::list2(...)
  verbose <- dots$verbose %||%
    SigBridgeRUtils::getFuncOption("verbose") %||%
    TRUE
  seed <- dots$seed %||% SigBridgeRUtils::getFuncOption("seed") %||% 123L
  assay <- dots$assay %||% "RNA"

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green("Start scSurvival screening"))
  }

  # * handling user parameters
  scsurvival_params <- scSurvivalParamSet(scsurvival_params = scsurvival_params)
  env_params <- scSurvivalEnvSet(env_params = env_params)

  python <- FindPy(
    env_params = env_params,
    method_name = "scSurvival",
    verbose = verbose
  )

  if (verbose) {
    ts_cli$cli_alert_info("Check Python execution environment")
  }
  reticulate::use_python(python)

  sc_expr <- SeuratObject::LayerData(sc_data, layer = "data", assay = assay)

  # Return result in expected format
  list(
    scRNA_data = modified_sc_data
  )
}

#' @keywords internal
#' @family scSurvival
scSurvivalEnvSet <- function(env_params = list()) {
  default <- list(
    env.name = "r-reticulate-scsurvival-nvidia",
    env.type = "conda",
    env.method = "environment",
    env.file = system.file(
      "conda/scSurvival_nvidia_environment.yml",
      package = "SigBridgeR"
    ),
    env.python_verion = "3.12.13",
    env.packages = c(
      "numpy" = "1.26.4"
      # and more
    ),
    env.recreate = FALSE,
    env.use_conda_forge = TRUE,
    env.verbose = getFuncOption("verbose")
  )

  utils::modifyList(default, env_params)
}

scSurvivalParamSet <- function(scsurvival_params = list()) {
  default <- list()
  utils::modifyList(default, scsurvival_params)
}
