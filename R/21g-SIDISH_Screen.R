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
#' @param sidish_params List of SIDISH algorithm parameters including:
#'   **Preprocessing parameters:**
#'   - `patient_id`: column name for patient identifier in metadata (default: `"Sample"`)
#'   - `celltype_name`: column name for cell type annotation in metadata (default: `"celltype_major"`)
#'   - `processed`: whether input data is already preprocessed (default: `TRUE`)
#'   - `n_genes_by_counts`: minimum number of genes expressed per cell filter threshold (default: `5000`)
#'   - `pct_counts_mt`: maximum percentage of mitochondrial genes filter threshold (default: `10`)
#'   - `batch_correction`: whether to perform batch correction (default: `FALSE`)
#'   - `survival_`: column name for survival time in phenotype data (default: `"time"`)
#'   - `status`: column name for event status in phenotype data (default: `"status"`)
#'
#'   **Execution environment:**
#'   - `device`: computation device, `"cuda"` for GPU acceleration or `"cpu"` for CPU-only (default: `"cuda"`)
#'   - `use_spatial_graph`: whether to use spatial graph information (default: `FALSE`)
#'   - `k_neighbors`: number of neighbors for graph construction (default: `NULL`, auto-detected)
#'
#'   **Phase 1: VAE training parameters:**
#'   - `phase1_epochs`: total epochs for VAE training (default: `225`)
#'   - `phase1_i_epochs`: interval epochs for VAE intermediate evaluation (default: `20`)
#'   - `phase1_latent_size`: dimensionality of latent space (default: `32`)
#'   - `phase1_layer_dims`: hidden layer dimensions as integer vector (default: `c(512, 128)`)
#'   - `phase1_batch_size`: batch size for VAE training (default: `256`)
#'   - `phase1_optimizer`: optimizer algorithm (default: `"Adam"`)
#'   - `phase1_lr`: learning rate for VAE encoder/decoder (default: `1e-4`)
#'   - `phase1_lr_3`: learning rate for additional VAE component (default: `1e-4`)
#'   - `phase1_dropout`: dropout rate for VAE layers (default: `0`)
#'   - `phase1_type`: VAE layer type, `"Dense"` or `"Normal"` (default: `"Dense"`)
#'
#'   **Phase 2: Deep Cox training parameters:**
#'   - `phase2_epochs`: total epochs for Cox model training (default: `500`)
#'   - `phase2_hidden`: number of hidden units in Cox model (default: `128`)
#'   - `phase2_lr`: learning rate for Cox model (default: `1e-4`)
#'   - `phase2_dropout`: dropout rate for Cox model (default: `0`)
#'   - `phase2_test_size`: proportion of data held out for testing (default: `0.2`)
#'   - `phase2_batch_size_bulk`: batch size for bulk data in Cox training (default: `256`)
#'
#'   **Training & risk definition parameters:**
#'   - `train_iterations`: number of risk score iteration rounds (default: `5`)
#'   - `train_percentile`: percentile threshold for high-risk cell selection (default: `0.95`)
#'   - `train_steepness`: steepness parameter for risk score transformation (default: `30`)
#'   - `train_path`: directory path for saving intermediate results (default: `"./SIDISH_res/"`)
#'   - `train_num_workers`: number of data loading workers (default: `0`)
#'   - `train_distribution_fit`: distribution fitting method, `"fitted"` or `"default"` (default: `"fitted"`)
#' @param ... Additional arguments passed to the function. Common parameters include:
#'   \describe{
#'     \item{verbose}{Logical. Whether to print verbose output (default: `TRUE`).}
#'     \item{seed}{Integer. Random seed for reproducibility (default: `123L`).}
#'     \item{assay}{Character. Assay to use for screening (default: `"RNA"`).}
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
  label_type = "SIDISH",
  phenotype_class = "survival",
  sidish_params = list(),
  ...
) {
  check_installed("rSIDISH", action = \(pkg, ...) {
    check_installed("pak")
    pak::pak("Exceret/rSIDISH")
  })

  chk::chk_is(sc_data, "Seurat")
  # Validate phenotype_class parameter
  if (phenotype_class != "survival") {
    cli::cli_abort(c("x" = "Currently phenotype_class must be 'survival'"))
  }

  # Extract additional arguments
  dots <- list(...)
  verbose <- dots$verbose %||% getFuncOption("verbose") %||% TRUE
  seed <- dots$seed %||% getFuncOption("seed") %||% 123L
  assay <- dots$assay %||% "RNA"
  env_params <- dots$env_params

  if (is.list(env_params)) {
    rlang::check_installed("lifecycle")
    lifecycle::deprecate_warn(
      "3.8.1",
      "DoSIDISH(env_params)",
      details = "After the R-side binding of SIDISH is updated to version >=0.0.2, user-provided environment is no longer required."
    )
  }

  # * handling user parameters
  sidish_params <- SIDISHParamSet(sidish_params = sidish_params)

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
    !!!sidish_params,
    seed = seed
  )
  # A Seurat
  res <- AddMisc(res, SIDISH_paramss = sidish_params, SIDISH_type = label_type)

  # Return result in expected format
  list(
    scRNA_data = res
  )
}

#' SIDISH Env Setting
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' @keywords internal
#' @family SIDISH
SIDISHEnvSet <- function(env_params = list(), device = c("cuda", "cpu")) {
  rlang::check_installed("lifecycle")
  lifecycle::deprecate_warn(
    "3.8.1",
    "SIDISHEnvSet()",
    details = "After the R-side binding of SIDISH is updated to version >=0.0.2, user-provided environment is no longer required. This function has been deprecated and is retained for reference only"
  )
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
SIDISHParamSet <- function(sidish_params = list()) {
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

  utils::modifyList(default, sidish_params)
}
