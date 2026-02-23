#' @title Perform SIDISH Screening Analysis
#'
#' @param matched_bulk Matrix or data frame of preprocessed bulk RNA-seq expression
#'        data (genes x samples). Column names must match names/IDs in `phenotype`.
#' @param sc_data A Seurat object containing scRNA-seq data to be screened.
#' @param phenotype Phenotype data, either:
#'        - Patient survival Data frame with row names matching `matched_bulk` columns, colnames named "time" and "status"
#' @param label_type Character specifying phenotype label type
#' @param phenotype_class Type of phenotypic outcome (must be consistent with input data):
#'        - `"survival"`: Survival infomation
#' @param assay Seurat assay name, default: `"RNA"`.
#' @param sidish_param Parameters adjusting SIDISH. Use `SigBridgeR:::SIDISHParamSet()` to see all available parameters and fetch default values.
#' @param env_params Parameters adjusting python environment. Use `SigBridgeR:::SIDISHEnvSet()`to see all available parameters and fetch default values.
#' @param ... Additional arguments passed to the function. Common parameters include:
#'   \describe{
#'     \item{verbose}{Logical. Whether to print verbose output (default: `TRUE`).}
#'   }
#'
#' @return A named list containing:
#'   \describe{
#'     \item{scRNA_data}{Modified single-cell data object with integrated screening results.}
#'   }
#' @family screen_method
#' @family SIDISH
#' @export
DoSIDISH <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = NULL,
  phenotype_class = "survival",
  assay = "RNA",
  sidish_param = list(),
  env_params = list(),
  ...
) {
  chk::chk_is(sc_data, "Seurat")
  # Validate phenotype_class parameter
  if (phenotype_class != "survival") {
    cli::cli_abort(c("x" = "Currently phenotype_class must be 'survival'"))
  }

  # Extract additional arguments
  dots <- list(...)
  verbose <- dots$verbose %||% TRUE
  label_type <- label_type %||% "SIDISH"

  # * handling user parameters
  sidish_param <- SIDISHParamSet(sidish_param = sidish_param)
  env_params <- SIDISHEnvSet(
    env_params = env_params,
    device = sidish_param$device
  )

  if (!rSIDISH::detect_gpu(verbose = FALSE) && sidish_param$device == "cuda") {
    cli::cli_abort(c(
      "x" = "GPU not detected, please set {.pkg sidish_param = list(device = 'cpu')}"
    ))
  }

  python <- FindPy(
    env_params = env_params,
    method_name = "SIDISH",
    verbose = verbose
  )

  reticulate::use_python(python)

  res <- rSIDISH::sidish(
    matched_bulk = matched_bulk,
    sc_data = sc_data,
    phenotype = phenotype,
    label_type = label_type,
    phenotype_class = phenotype_class,
    assay = assay,
    verbose = verbose,
    sidish_tools = system.file(
      "python/01_training_SIDISH.py",
      package = "rSIDISH"
    ),
    !!!sidish_param,
    seed = seed
  )
  # A Seurat
  res <- AddMisc(res, SIDISH_params = sidish_param)

  # Return result in expected format
  list(
    scRNA_data = res
  )
}

#' @keywords internal
#' @family SIDISH
SIDISHEnvSet <- function(env_params = list(), device = c("cuda", "cpu")) {
  default <- list(
    env.name = if (device == "cuda") {
      "r-reticulate-sidish-nvidia"
    } else {
      "r-reticulate-sidish-cpu"
    },
    env.type = "conda",
    env.method = "environment",
    env.file = if (device == "cuda") {
      system.file(
        "conda/SIDISH_nvidia_environment.yml",
        package = "SigBridgeR"
      )
    } else {
      system.file(
        "conda/SIDISH_cpu_environment.yml",
        package = "SigBridgeR"
      )
    },
    env.python_verion = "3.12.12",
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

#' @keywords internal
#' @family SIDISH
SIDISHParamSet <- function(sidish_param = list()) {
  default <- list(
    # Preprocessing parameters
    patient_id = "Sample",
    celltype_name = "celltype_major",
    processed = TRUE,
    n_genes_by_counts = 5000L,
    pct_counts_mt = 10L,
    batch_correction = FALSE,
    survival_ = "time",
    status = "status",

    # Execution environment
    device = "cuda", # or "cpu"
    use_spatial_graph = FALSE,
    k_neighbors = NULL,

    # Phase 1: VAE training
    phase1_epochs = 225L,
    phase1_i_epochs = 20L,
    phase1_latent_size = 32L,
    phase1_layer_dims = c(512L, 128L), # R vector -> py list
    phase1_batch_size = 256L,
    phase1_optimizer = "Adam",
    phase1_lr = 1e-4,
    phase1_lr_3 = 1e-4,
    phase1_dropout = 0L,
    phase1_type = "Dense", # or "Normal"

    # Phase 2: Deep Cox training
    phase2_epochs = 500L,
    phase2_hidden = 128L,
    phase2_lr = 1e-4,
    phase2_dropout = 0L,
    phase2_test_size = 0.2,
    phase2_batch_size_bulk = 256L,

    # Training & risk definition
    train_iterations = 5L,
    train_percentile = 0.95,
    train_steepness = 30L,
    train_path = "./SIDISH_res/",
    train_num_workers = 0L,
    train_distribution_fit = "fitted" # or "default"
  )

  utils::modifyList(default, sidish_param)
}
