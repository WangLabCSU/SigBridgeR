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
#' @param tmp_dir Temporary directory for intermediate files (default: `"tmp"`)
#' @param env_params List of environment parameters for Python setup including:
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
  tmp_dir = "tmp",
  # DEGAS environment
  env_params = list(),
  # DEGAS parameters
  degas_params = list(),
  normality_test_method = c(
    "jarque-bera",
    "d'agostino",
    "kolmogorov-smirnov"
  ),
  ...
) {
  CheckInstalled("Exceret/DEGAS")
  # * inputs check
  chk::chk_range(select_fraction)
  chk::chk_is(matched_bulk, c("matrix", "data.frame"))
  chk::chk_not_any_na(matched_bulk)
  if (!chk::vld_is(sc_data, c("Seurat", "Matrix", "matrix"))) {
    cli::cli_abort(c(
      "x" = "{.arg sc_data} cannot be of type {.cls {class(sc_data)}}",
      ">" = "Available types: {.cls {c('Seurat', 'Matrix', 'matrix')}}"
    ))
  }
  chk::chk_character(label_type)
  if (!is.null(sc_data.pheno_colname)) {
    chk::chk_is(sc_data.pheno_colname, c("character"))
  }
  chk::chk_list(env_params)
  chk::chk_list(degas_params)
  phenotype_class <- SigBridgeRUtils::MatchArg(
    phenotype_class,
    c("binary", "continuous", "survival"),
    NULL
  )

  dots <- rlang::list2(...)
  verbose <- dots$verbose %||%
    SigBridgeRUtils::getFuncOption("verbose") %||%
    TRUE
  assay <- dots$assay %||% "RNA"

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green("Starting DEGAS Screen"))
  }

  # Auto-choose normality test method
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
    cli::cli_abort(c(
      "x" = "{.arg tmp_dir} must be specified."
    ))
  }
  if (!grepl("/$", tmp_dir)) {
    tmp_dir <- paste0(tmp_dir, "/")
  }

  env_params <- DEGASEnvSet(user_list = env_params)
  degas_params <- DEGASParamSet(user_list = degas_params)

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

  # formatting phenotype
  pheno_df <- if (!is.null(phenotype)) {
    switch(
      phenotype_class,
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
  sc_mat <- if (inherits(sc_data, "Seurat")) {
    SeuratObject::LayerData(sc_data, layer = "data", assay = assay)
  } else {
    sc_data
  }
  cm_genes <- intersect(rownames(matched_bulk), rownames(sc_mat))

  if (length(cm_genes) == 0) {
    cli::cli_abort(c(
      "x" = "No common genes found between single cell data and bulk data",
      ">" = "Please check the inputs"
    ))
  }

  t_sc_mat <- Matrix::t(sc_mat[cm_genes, ])
  t_matched_bulk <- Matrix::t(matched_bulk[cm_genes, ])

  if (verbose) {
    ts_cli$cli_alert_info("Setting up Environment...")
  }

  # * Check if environment exists and locate python
  if (is.null(degas_params$DEGAS.pyloc)) {
    degas_params$DEGAS.pyloc <- FindPy(
      env_params = env_params,
      verbose = verbose,
      !!!dots
    )
  }

  rlang::check_installed("reticulate")
  reticulate::use_python(
    degas_params$DEGAS.pyloc,
    required = TRUE
  )

  if (verbose) {
    ts_cli$cli_alert_info("Training DEGAS model...")
  }

  # DEGAS needs some global variables to be set up
  list2env(degas_params, envir = .GlobalEnv)
  on.exit({
    # Clean up global variables
    rm(list = names(degas_params), envir = .GlobalEnv)
    gc(verbose = FALSE)
  })

  # Train DEGAS model
  ccModel1 <- DEGAS::runCCMTLBag.optimized(
    scExp = t_sc_mat,
    scLab = sc_pheno,
    patExp = t_matched_bulk,
    patLab = pheno_df,
    tmpDir = tmp_dir,
    model_type = DEGAS.model_type,
    architecture = DEGAS.architecture,
    FFdepth = DEGAS.ff_depth,
    Bagdepth = DEGAS.bag_depth,
    DEGAS.seed = DEGAS.seed,
    verbose = verbose
  )
  names(ccModel1) <- glue::glue("ccModel_{seq_len(DEGAS.bag_depth)}")

  if (verbose) {
    ts_cli$cli_alert_info("Predicting and Labeling")
  }

  # Predict with DEGAS model
  t_sc_preds <- data.table::as.data.table(DEGAS::predClassBag.optimized(
    ccModel = ccModel1,
    Exp = t_sc_mat,
    scORpat = 'pat'
  ))

  if (phenotype_class == "survival") {
    data.table::setnames(t_sc_preds, "Hazard")
  } else {
    pheno_df_colnames <- if (!is.null(pheno_df)) {
      colnames(pheno_df)
    } else {
      NULL
    }
    if (!is.null(pheno_df_colnames)) {
      data.table::setnames(t_sc_preds, pheno_df_colnames)
    }
  }
  t_sc_preds[, "cell_id" := rownames(t_sc_mat)]

  if (verbose) {
    ts_cli$cli_alert_info("Labeling screened cells")
  }

  # What we get from DEGAS are the probabilities of the cells belonging to each class.
  # We need to convert these to labels.
  t_sc_preds <- switch(
    phenotype_class,
    "binary" = DEGAS::LabelBinaryCells(
      pred_dt = t_sc_preds,
      pheno_colnames = pheno_df_colnames,
      select_fraction = select_fraction,
      test_method = normality_test_method,
      min_threshold = min_thresh,
      verbose = verbose
    ),
    "continuous" = DEGAS::LabelContinuousCells(
      pred_dt = t_sc_preds,
      verbose = verbose
    ),
    "survival" = DEGAS::LabelSurvivalCells(
      pred_dt = t_sc_preds,
      select_fraction = select_fraction,
      test_method = normality_test_method,
      min_threshold = min_thresh,
      verbose = verbose
    )
  )

  meta2add <- t_sc_preds[,
    !names(t_sc_preds) %in% c("cell_id", "diff"),
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
  ) %>%
    SigBridgeRUtils::AddMisc(
      DEGAS_type = label_type,
      DEGAS_para = degas_params,
      cover = FALSE
    )
  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green("DEGAS Screen done."))
  }

  # result
  list(
    scRNA_data = sc_data,
    model = ccModel1,
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
  call = rlang::caller_env()
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
    cli::cli_abort(
      c(
        "x" = "Please specify at least one phenotype, currently both are {.val NULL}"
      ),
      class = "NoPhenotypeProvided",
      call = call
    )
  }
  paste0(model_type.first, model_type.last)
}

#' @title Set Default Values for Environment Configuration.
#'
#' @keywords internal
#' @family DEGAS
DEGASEnvSet <- function(user_list) {
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

#' @title Find Python Location from Environment Setting
#' @description
#' Locate the Python executable path from a predefined environment configuration list;
#' if no suitable environment exists, create a new Python environment and return its path.
#'
#' @keywords internal
#'
FindPy <- function(env_params, method_name = "DEGAS", verbose = TRUE, ...) {
  # * Check if environment exists
  existing_envs <- ListPyEnv(
    env_type = env_params$env.type,
    verbose = FALSE
  )
  if (
    !env_params$env.name %chin% existing_envs$name ||
      env_params$env.recreate
  ) {
    choice <- utils::askYesNo(paste0(
      "Create a new conda environment for ",
      method_name,
      "?"
    ))

    if (!isTRUE(choice)) {
      cli::cli_abort(c(
        "x" = "Aborted."
      ))
    }

    # setup environment or recreate if needed
    rlang::exec(
      SetupPyEnv,
      env_type = env_params$env.type,
      env_name = env_params$env.name,
      python_version = env_params$env.python_verion,
      packages = env_params$env.packages,
      recreate = env_params$env.recreate,
      verbose = env_params$env.verbose,
      !!!if (env_params$env.type == "conda") {
        list(
          use_conda_forge = env_params$env.use_conda_forge,
          env_file = env_params$env.file,
          method = env_params$env.method
        )
      },
      ...
    )
  } else if (verbose) {
    ts_cli$cli_alert_info(
      "Existing environment {.val {env_params$env.name}} found"
    )
  }

  existing_envs$python[existing_envs$name == env_params$env.name]
}
