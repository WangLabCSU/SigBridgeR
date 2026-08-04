# ---- 2. Do DEGAS ----

#' Validate and Prepare DoDEGAS Parameters
#'
#' @description
#' Internal helper that handles package installation checks, input validation,
#' default-value resolution, and cache configuration for [DoDEGAS()].
#'
#' @param select_fraction,min_thresh,matched_bulk,sc_data,phenotype
#'   Forwarded from [DoDEGAS()].
#' @param sc_data.pheno_colname,label_type,phenotype_class,tmp_dir
#'   Forwarded from [DoDEGAS()].
#' @param env_params,degas_params,normality_test_method
#'   Forwarded from [DoDEGAS()].
#' @param ... Additional dots forwarded from [DoDEGAS()].
#'
#' @return A named list with elements: `phenotype_class`, `normality_test_method`,
#'   `verbose`, `assay`, `load_cache`, `save_cache`, `env_params`, `degas_params`,
#'   `cache_config`, `dots`.
#'
#' @keywords internal
#' @family DEGAS
ValidateDEGASParams <- function(
  select_fraction,
  min_thresh,
  matched_bulk,
  sc_data,
  phenotype,
  sc_data.pheno_colname,
  label_type,
  phenotype_class,
  tmp_dir,
  env_params = lifecycle::deprecated(),
  degas_params,
  normality_test_method,
  ...
) {
  # TODO Update DEGAS implemention
  # -- package checks -------------------------------------------------------
  check_installed("DEGAS", action = \(pkg, ...) {
    check_installed("pak")
    pak::pak("Exceret/DEGAS")
  })

  # -- input validation -----------------------------------------------------
  chk::chk_range(select_fraction)
  chk::chk_is(matched_bulk, c("matrix", "data.frame"))
  chk::chk_not_any_na(matched_bulk)
  if (!chk::vld_is(sc_data, c("Seurat", "Matrix", "matrix"))) {
    Abort(
      "{.arg sc_data} cannot be of type {.cls {class(sc_data)}}",
      "Available types: {.cls {c('Seurat', 'Matrix', 'matrix')}}"
    )
  }
  chk::chk_character(label_type)
  if (!is.null(sc_data.pheno_colname)) {
    chk::chk_is(sc_data.pheno_colname, c("character"))
  }
  chk::chk_list(degas_params)
  phenotype_class <- SigBridgeRUtils::MatchArg(
    phenotype_class,
    c("binary", "continuous", "survival"),
    NULL
  )

  normality_test_method <- SigBridgeRUtils::MatchArg(
    normality_test_method,
    c(
      "jarque-bera",
      "d'agostino",
      "kolmogorov-smirnov"
    )
  )

  # DEGAS path must contain "/"
  if (is.null(tmp_dir)) {
    Abort(
      "{.arg tmp_dir} must be specified."
    )
  }
  if (!endsWith(tmp_dir, "/")) {
    tmp_dir <- paste0(tmp_dir, "/")
  }

  # -- resolve env & degas defaults -----------------------------------------
  degas_params <- DEGASParamSet(user_list = degas_params)
  if (lifecycle::is_present(env_params)) {
    lifecycle::deprecate_warn(
      "4.0.0",
      "DoDEGAS(env_params = )"
    )
  }

  degas_params$DEGAS.model_type <- if (
    length(degas_params$DEGAS.model_type) != 1
  ) {
    # model.type auto-detection
    DEGASModelDetect(
      sc_pheno_colname = sc_data.pheno_colname,
      bulk_pheno = phenotype
    )
  } else {
    # check validity
    SigBridgeRUtils::MatchArg(
      degas_params$DEGAS.model_type,
      c('BlankClass', 'ClassBlank', 'ClassClass', 'ClassCox', 'BlankCox'),
      NULL
    )
  }

  # Architecture auto-chosen
  degas_params$DEGAS.architecture <- SigBridgeRUtils::MatchArg(
    degas_params$DEGAS.architecture,
    c(
      "DenseNet", # default
      "Standard"
    )
  )

  # -- process dots ---------------------------------------------------------
  dots <- list2(...)
  verbose <- dots$verbose %||%
    SigBridgeRUtils::getFuncOption("verbose") %||%
    TRUE
  assay <- dots$assay %||% "RNA"
  load_cache <- dots$load_cache
  save_cache <- dots$save_cache %||% tmp_dir

  sc_mat <- if (inherits(sc_data, "Seurat")) {
    SeuratObject::LayerData(sc_data, layer = "data", assay = assay)
  } else {
    sc_data
  }
  cm_genes <- intersect(rownames(matched_bulk), rownames(sc_mat))

  if (length(cm_genes) == 0) {
    Abort(
      "No common genes found between single cell data and bulk data",
      "Please check the inputs"
    )
  }

  t_sc_mat <- Matrix::t(sc_mat[cm_genes, ])

  # -- build cache config ---------------------------------------------------
  cache_config <- ScreenMethodConfig(
    method_name = "DEGAS",
    param = get_env_vars(exclude = c("matched_bulk", "sc_data", "phenotype")),
    phenotype_class = phenotype_class,
    label_type = label_type
  )

  get_env_vars(exclude = c("matched_bulk", "sc_data", "phenotype"))
}

#' Train DEGAS Model and Optionally Cache
#'
#' @description
#' Internal helper that sets up the Python environment, trains the DEGAS model
#' using bootstrap-aggregated multi-task learning, and optionally persists the
#' result to cache.
#'
#' @param matched_bulk,sc_data,phenotype,sc_data.pheno_colname,label_type
#'   Forwarded from [DoDEGAS()].
#' @param p Named list of validated parameters from [ValidateDEGASParams()].
#'
#' @return The trained `ccModel1` object (list of DEGAS models).
#'
#' @keywords internal
#' @family DEGAS
TrainDEGASModel <- function(
  matched_bulk,
  sc_data,
  phenotype,
  t_sc_mat,
  cm_genes,
  sc_data.pheno_colname,
  label_type,
  p
) {
  if (p$verbose) {
    ts_cli$cli_alert_info("Setting up Environment")
  }

  # * Check if environment exists and locate python
  check_installed("reticulate")
  p$degas_params$DEGAS.pyloc <- if (
    !reticulate::py_available(initialize = FALSE)
  ) {
    Abort(
      "Python is not available.",
      "Please configure python executable path with `reticulate::use_python(path)`"
    )
  } else {
    py <- reticulate::py_config()$python
    if (p$verbose) {
      ts_cli$cli_alert_info("Found {.file {.path {py}}")
    }
    py
  }

  # -- train DEGAS model --------------------------------------------------
  if (p$verbose) {
    ts_cli$cli_alert_info("Training DEGAS model")
  }

  # DEGAS needs some global variables to be set up
  on.exit({
    # Clean up global variables
    gc(verbose = FALSE)
  })

  # formatting phenotype
  pheno_df <- if (!is.null(phenotype)) {
    switch(
      p$phenotype_class,
      "binary" = DEGAS::Vec2sparse(
        phenotype,
        col_prefix = label_type
      ),
      "continuous" = DEGAS::Vec2sparse(
        phenotype,
        col_prefix = label_type
      ),
      "survival" = as.matrix(phenotype),
    )
  } else {
    NULL
  }

  # Check if single-cell level phenotype is specified
  meta_data <- sc_data[[]]
  sc_pheno <-
    if (
      !is.null(sc_data.pheno_colname) &&
        sc_data.pheno_colname %chin% colnames(meta_data)
    ) {
      DEGAS::Vec2sparse(meta_data[[sc_data.pheno_colname]])
    } else {
      if (!is.null(sc_data.pheno_colname)) {
        cli::cli_warn(
          "Single-cell phenotype specified but not found in meta.data, using {.val NULL}"
        )
      }
      NULL
    }

  # anndata-like data formats
  t_matched_bulk <- Matrix::t(matched_bulk[cm_genes, ])

  # Train DEGAS model
  ccModel1 <- DEGAS::runCCMTLBag.optimized(
    scExp = t_sc_mat,
    scLab = sc_pheno,
    patExp = t_matched_bulk,
    patLab = pheno_df,
    tmpDir = p$tmp_dir,
    model_type = p$degas_paramsDEGAS.model_type,
    architecture = p$degas_paramsDEGAS.architecture,
    FFdepth = p$degas_paramsDEGAS.ff_depth,
    Bagdepth = p$degas_paramsDEGAS.bag_depth,
    DEGAS.pyloc = p$degas_params$DEGAS.pyloc,
    DEGAS.toolsPath = p$degas_params$DEGAS.toolsPath,
    DEGAS.train_steps = p$degas_params$DEGAS.train_steps,
    DEGAS.scbatch_sz = p$degas_params$DEGAS.scbatch_sz,
    DEGAS.patbatch_sz = p$degas_params$DEGAS.patbatch_sz,
    DEGAS.hidden_feats = p$degas_params$DEGAS.hidden_feats,
    DEGAS.do_prc = p$degas_params$DEGAS.do_prc,
    DEGAS.lambda1 = p$degas_params$DEGAS.lambda1,
    DEGAS.lambda2 = p$degas_params$DEGAS.lambda2,
    DEGAS.lambda3 = p$degas_params$DEGAS.lambda3,
    DEGAS.seed = p$degas_params$DEGAS.seed,
    verbose = p$verbose,
    force_rewrite = p$force_rewrite %||% FALSE,
  )
  names(ccModel1) <- paste0(
    "ccModel_",
    seq_len(p$degas_params$DEGAS.bag_depth)
  )

  # -- save mode: persist ccModel1 to cache --------------------------------
  if (!is.null(p$save_cache)) {
    cache <- ScreenMethodCache(
      cache_path = p$save_cache,
      cache_config_path = file.path(p$save_cache, "cache_config.json"),
      cache_data = list(ccModel1 = ccModel1, pheno_df = pheno_df),
      screen_method_config = p$cache_config
    )
    CacheSysCall(
      mode = "save",
      path = p$load_cache,
      cache = cache,
      verbose = p$verbose,
      timestamp = p$dots$timestamp
    )
  }
  list(ccModel1 = ccModel1, pheno_df = pheno_df)
}

#' @title Run DEGAS Analysis for Single-Cell and Bulk RNA-seq Data Integration
#'
#' @description
#' This function performs DEGAS to integrate
#' single-cell and bulk RNA-seq data, identifying phenotype-associated cells using
#' a bootstrap aggregated multi-task learning approach.
#'
#' @param select_fraction The top percentage of selected cells will be considered as Positive cells, without considering how much larger the possible correlation coefficient of the observation group is compared to that of the control group. Only usedl when `phenotype_class` is "binary" or "survival". (default: `0.05`)
#' @param min_thresh DEGAS will calculate the possible correlation coefficients for each cell related to the phenotype. When the coefficient of the observation group is at least `min_thresh` larger than that of the control group, it can be considered related to the phenotype and will be marked as Positive. The priority of `min_thresh` is higher than that of `select_fraction.` (default: `0.4`)
#' @param matched_bulk Bulk RNA-seq data as matrix or data.frame (rows=genes, columns=samples)
#' @param sc_data Single-cell data as Seurat object containing RNA assay
#' @param phenotype Bulk-level phenotype data. For classification: binary matrix with one-hot encoding.
#'   For survival: matrix with two columns (time and event status). Can be NULL, matrix, data.frame, or vector.
#' @param sc_data.pheno_colname Column name for single-cell phenotype in metadata (if available), default: `NULL`
#' @param label_type Label type for DEGAS results (default: "`DEGAS"`)
#' @param phenotype_class Type of phenotype: "binary" (classification), "continuous", or "survival"
#' @param tmp_dir (Soft-deprecated) Temporary directory for DEGAS intermediate
#'   files (default: `"DEGAS_res"`). Acts as fallback for \code{load_cache} and
#'   \code{save_cache} when not specified via \code{...}. Prefer using
#'   \code{load_cache} / \code{save_cache} in \code{...} instead. See [CacheSetHere()].
#' @param env_params `r lifecycle::badge("deprecated")` List of environment parameters for Python setup including:
#'   - env.name: environment name (default: `"r-reticulate-degas"`)
#'   - env.type: environment type "conda", "environment", or "venv" (default: `"environment"`)
#'   - env.method: environment setup method "system", "conda" (default: `"system"`)
#'   - env.file: path to environment file (default: system.file("conda/DEGAS_environment.yml", package = "SigBridgeR"))
#'   - env.python_version: Python version (default: "3.9.15")
#'   - env.packages: named vector of Python packages and versions (default: c("tensorflow" = "2.4.1", "protobuf" = "3.20" ,"numpy" = "any"))
#'   - env.recreate: whether to recreate environment (default: `FALSE`)
#'   - env.use_conda_forge: whether to use conda-forge channel (conda only, default: `TRUE`)
#'   - env.verbose: verbose output (default: `FALSE`)
#' @param degas_params List of DEGAS algorithm parameters including:
#'   - DEGAS.model_type: model type ("BlankClass", "ClassBlank", "ClassClass", "ClassCox", "BlankCox")
#'   - DEGAS.architecture: "Standard" (feed forward) or "DenseNet" (dense net), default: `"DenseNet"`
#'   - DEGAS.ff_depth: number of layers in model (>=1, default: `3`)
#'   - DEGAS.pyloc: path to Python executable (default: `NULL`, automatic detection)
#'   - DEGAS.bag_depth: bootstrap aggregation depth (>=1, default: `5`)
#'   - DEGAS.train_steps: training steps (default: `2000`)
#'   - DEGAS.scbatch_sz: single-cell batch size (default: `200`)
#'   - DEGAS.patbatch_sz: patient batch size (default: `50`)
#'   - DEGAS.hidden_feats: hidden features (default: `50`)
#'   - DEGAS.do_prc: dropout percentage (default: `0.5`)
#'   - DEGAS.lambda1: regularization parameter 1 (default: `3.0`)
#'   - DEGAS.lambda2: regularization parameter 2 (default: `3.0`)
#'   - DEGAS.lambda3: regularization parameter 3 (default: `3.0`)
#'   - DEGAS.seed: random seed (default: `2`)
#' @param normality_test_method Method for normality testing: "jarque-bera", "d'agostino", or "kolmogorov-smirnov"
#' @param ... Additional arguments. Currently supports:
#'    - `verbose`: Logical indicating whether to print progress messages. Defaults to `TRUE`.
#'    - `assay`: Name of assay to use. Defaults to "RNA".
#'    - `load_cache`: Cache directory path for loading a precomputed DEGAS model.
#'      Supports root-level, cache-level, or parent-level paths. See [CacheSetHere()].
#'    - `save_cache`: Cache directory path for saving the trained DEGAS model.
#'      Supports root-level or parent-level paths. See [CacheSetHere()].
#'
#' @return A list containing:
#'   - scRNA_data: Seurat object with DEGAS labels added to metadata
#'   - model: The model trained using the input data, andit can be used for cell classification prediction.
#'   - DEGAS_prediction: Data table with DEGAS predictions containing:
#'     * Predicted label probabilities for each cell
#'     * Cell labels ("Positive"/"Other") based on selection criteria
#'     * Difference scores for binary phenotypes
#'     * Cell identifiers
#'
#' @details
#' The function performs the following steps:
#' 1. Validates input data and parameters
#' 2. Sets up Python environment with required dependencies
#' 3. Trains bootstrap aggregated DEGAS model using `runCCMTLBag`
#' 4. Generates cell-level predictions using `predClassBag`
#' 5. Applies statistical testing to identify phenotype-associated cells
#' 6. Labels cells as "Positive" or "Other" based on selection criteria
#'
#' Model type is automatically determined:
#' - BlankClass: only bulk phenotype specified (scLab = NULL)
#' - ClassBlank: only single-cell phenotype specified (patLab = NULL)
#' - ClassClass: both single-cell and bulk phenotypes specified
#' - ClassCox: single-cell phenotype + bulk survival data
#' - BlankCox: only bulk survival data specified
#'
#' @examples
#' \dontrun{
#' # Binary classification example
#' result <- DoDEGAS(
#'   select_fraction = 0.05, # `select_fraction` only used in binary and survival phenotyping
#'   matched_bulk = bulk_matrix,
#'   sc_data = seurat_obj,
#'   phenotype = bulk_phenotype,
#'   phenotype_class = "binary"
#' )
#'
#' # Survival analysis example
#' result <- DoDEGAS(
#'   select_fraction = 0.05, # `select_fraction` only used in binary and survival phenotyping
#'   matched_bulk = bulk_matrix,
#'   sc_data = seurat_obj,
#'   phenotype = survival_data,
#'   phenotype_class = "survival"
#' )
#' }
#'
#' @export
#' @family screen_method
#' @family DEGAS
#'
#' @references Johnson TS, Yu CY, Huang Z, Xu S, Wang T, Dong C, et al. Diagnostic Evidence GAuge of Single cells (DEGAS): a flexible deep transfer learning framework for prioritizing cells in relation to disease. Genome Med. 2022 Feb 1;14(1):11.
#'
DoDEGAS <- function(
  select_fraction = 0.05,
  min_thresh = 0.4,
  matched_bulk,
  sc_data,
  phenotype = NULL,
  sc_data.pheno_colname = NULL,
  label_type = "DEGAS",
  phenotype_class = c("binary", "continuous", "survival"),
  # A directory for intermediate files
  tmp_dir = "DEGAS_res",
  # DEGAS environment
  env_params = lifecycle::deprecated(),
  # DEGAS parameters
  degas_params = list(),
  normality_test_method = c(
    "jarque-bera",
    "d'agostino",
    "kolmogorov-smirnov"
  ),
  ...
) {
  # -- validate & prepare all parameters -----------------------------------
  p <- do.call(ValidateDEGASParams, c(get_env_vars(), list(...)))

  if (p$verbose) {
    ts_cli$cli_alert_info(cli::col_green("Starting DEGAS Screen"))
  }

  # -- load mode: restore cached ccModel1 ------------------------------------
  res <- if (!is.null(p$load_cache)) {
    cache <- CacheSysCall(
      mode = "load",
      path = p$load_cache,
      cache = p$cache_config,
      verbose = p$verbose,
      timestamp = p$dots$timestamp,
    )

    if (p$verbose) {
      ts_cli$cli_alert_info(
        cli::col_green("Loaded DEGAS cache from {.path {cache_dir}}")
      )
    }
    cache
  } else {
    TrainDEGASModel(
      matched_bulk = matched_bulk,
      sc_data = sc_data,
      phenotype = phenotype,
      t_sc_mat = p$t_sc_mat,
      cm_genes = p$cm_genes,
      sc_data.pheno_colname = sc_data.pheno_colname,
      label_type = label_type,
      p = p
    )
  } # !is.null(p$load_cache)

  if (p$verbose) {
    ts_cli$cli_alert_info("Predicting and labeling")
  }

  # Predict with DEGAS model
  t_sc_preds <- data.table::as.data.table(DEGAS::predClassBag.optimized(
    ccModel = res$ccModel1,
    Exp = p$t_sc_mat,
    scORpat = 'pat'
  ))

  if (p$phenotype_class == "survival") {
    data.table::setnames(t_sc_preds, "Hazard")
  } else {
    pheno_df_colnames <- if (!is.null(res$pheno_df)) {
      colnames(res$pheno_df)
    } else {
      NULL
    }
    if (!is.null(pheno_df_colnames)) {
      data.table::setnames(t_sc_preds, pheno_df_colnames)
    }
  }
  t_sc_preds[, "cell_id" := rownames(p$t_sc_mat)]

  if (p$verbose) {
    ts_cli$cli_alert_info("Labeling screened cells")
  }

  # What we get from DEGAS are the probabilities of the cells belonging to each class.
  # We need to convert these to labels.
  t_sc_preds <- switch(
    p$phenotype_class,
    "binary" = DEGAS::LabelBinaryCells(
      pred_dt = t_sc_preds,
      pheno_colnames = pheno_df_colnames,
      select_fraction = select_fraction,
      test_method = p$normality_test_method,
      min_threshold = min_thresh,
      verbose = p$verbose
    ),
    "continuous" = DEGAS::LabelContinuousCells(
      pred_dt = t_sc_preds,
      verbose = p$verbose
    ),
    "survival" = DEGAS::LabelSurvivalCells(
      pred_dt = t_sc_preds,
      select_fraction = select_fraction,
      test_method = p$normality_test_method,
      min_threshold = min_thresh,
      verbose = p$verbose
    )
  )

  meta2add <- t_sc_preds[,
    !names(t_sc_preds) %chin% c("cell_id", "diff"),
    with = FALSE
  ]
  data.table::setnames(
    meta2add,
    names(meta2add),
    data.table::fifelse(
      names(meta2add) == "label",
      "DEGAS",
      paste0("DEGAS_", names(meta2add))
    )
  )

  # Record screening results
  sc_data <- Seurat::AddMetaData(
    object = sc_data,
    metadata = meta2add
  )
  sc_data <- SigBridgeRUtils::AddMisc(
    seurat_obj = sc_data,
    DEGAS_type = label_type,
    DEGAS_para = p$degas_params,
    cover = FALSE
  )
  if (p$verbose) {
    ts_cli$cli_alert_info(cli::col_green("DEGAS Screen done."))
  }

  # result
  list(
    scRNA_data = sc_data,
    model = res$ccModel1,
    DEGAS_prediction = t_sc_preds
  )
}

#' @title Auto-select Model for DEGAS
#' @description
#' Automatically select the DEGAS model based on single-cell labels and bulk-level phenotype information.
#'
#' @keywords internal
#' @family DEGAS
DEGASModelDetect <- function(
  sc_pheno_colname,
  bulk_pheno,
  call = caller_env()
) {
  # model.type auto-detection
  model_type.first <- ifelse(
    !is.null(sc_pheno_colname),
    "Class",
    "Blank"
  )
  model_type.last <- if (is.null(bulk_pheno)) {
    "Blank"
  } else {
    ifelse(is.vector(bulk_pheno), "Class", "Cox")
  }
  if (model_type.first == "Blank" && model_type.last == "Blank") {
    Abort(
      "Please specify at least one phenotype, currently both are {.val NULL}"
    )
  }
  paste0(model_type.first, model_type.last)
}

#' @title Set Default Values for Environment Configuration.
#' @param user_list Configuration parameters provided by users. Override the default values.
#' @keywords internal
#' @family DEGAS
DEGASEnvSet <- function(user_list = list()) {
  lifecycle::deprecate_warn(
    "4.0.0",
    "DEGASEnvSet()",
    details = "This function provides a default environment configuration for DEGAS.\
    In SigBridgeR 3..y.z, the env is setted by default. Now, change to user-defined.\
    Left here for reference."
  )

  # default environment and DEGAS parameters
  default_env_params <- list(
    env.name = "r-reticulate-degas",
    env.type = "conda",
    env.method = "environment",
    env.file = system.file(
      "conda/DEGAS_environment.yml",
      package = "SigBridgeR"
    ),
    env.python_verion = "3.9.15",
    env.packages = c(
      "tensorflow" = "2.4.1",
      "protobuf" = "3.20.3"
    ),
    env.recreate = FALSE,
    env.use_conda_forge = TRUE,
    env.verbose = getFuncOption("verbose")
  )
  utils::modifyList(default_env_params, user_list)
}

#' @title Set Default Parameter Values for DEGAS
#'
#' @param user_list Configuration parameters provided by users. Override the default values.
#'
#' @keywords internal
#' @family DEGAS
DEGASParamSet <- function(user_list) {
  default_degas_params <- list(
    DEGAS.model_type = c(
      "BlankClass", # only bulk level phenotype specified
      "ClassBlank", # only single cell level phenotype specified
      "ClassClass", # when both single cell level phenotype and bulk level phenotype specified
      "ClassCox", # when both single cell level phenotype and bulk level survival data specified
      "BlankCox" # only bulk level survival data specified
    ),
    DEGAS.architecture = c(
      "DenseNet", # a dense net network
      "Standard" # a feed forward network
    ),
    DEGAS.ff_depth = 3L,
    DEGAS.bag_depth = 5L,
    path.data = '',
    path.result = '',
    DEGAS.pyloc = NULL, # location of python executable
    DEGAS.toolsPath = file.path(.libPaths()[1], "DEGAS/DEGAS_tools/"),
    DEGAS.train_steps = 2000L,
    DEGAS.scbatch_sz = 200L,
    DEGAS.patbatch_sz = 50L,
    DEGAS.hidden_feats = 50L,
    DEGAS.do_prc = 0.5,
    DEGAS.lambda1 = 3.0,
    DEGAS.lambda2 = 3.0,
    DEGAS.lambda3 = 3.0,
    DEGAS.seed = SigBridgeRUtils::getFuncOption("seed")
  )
  utils::modifyList(default_degas_params, user_list)
}
