# ? ---- 2. DO SCISSOR ----

#' Validate and Prepare DoScissor Parameters
#'
#' @description
#' Internal helper that handles deprecated argument warnings, input validation,
#' dots processing, and default-value resolution for [DoScissor()].
#'
#' @param matched_bulk,sc_data,phenotype,label_type,alpha,cutoff,family
#'   Forwarded from [DoScissor()].
#' @param reliability_test,cell_evaluation Forwarded from [DoScissor()].
#' @param path2load_scissor_cache,path2save_scissor_inputs Deprecated arguments.
#' @param ... Additional dots forwarded from [DoScissor()].
#'
#' @return A named list with elements: `family`, `verbose`, `seed`, `assay`,
#'   `load_cache`, `save_cache`, `reliability_test`, `cell_evaluation`,
#'   `label_type`, `phenotype_class_cache`, `dots`.
#'
#' @keywords internal
#' @family scissor
ValidateScissorParams <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type,
  phenotype_class = c("binary", "continuous", "survival"),
  alpha,
  cutoff,
  family,
  reliability_test,
  cell_evaluation,
  path2load_scissor_cache,
  path2save_scissor_inputs,
  ...
) {
  # -- deprecated-argument warnings -----------------------------------------
  if (lifecycle::is_present(path2load_scissor_cache)) {
    lifecycle::deprecate_warn(
      "4.0.0",
      "DoScissor(path2load_scissor_cache = )",
      "DoScissor(load_cache = )"
    )
  }
  if (lifecycle::is_present(path2save_scissor_inputs)) {
    lifecycle::deprecate_warn(
      "4.0.0",
      "DoScissor(path2save_scissor_inputs = )",
      "DoScissor(save_cache = )"
    )
  }
  if (lifecycle::is_present(family)) {
    lifecycle::deprecate_warn(
      "4.0.0",
      "DoScissor(family = )",
      "DoScissor(phenotype_class = )"
    )
    phenotype_class <- switch(
      family,
      "binomial" = "binary",
      "cox" = "survival",
      "gaussian" = "continuous",
      Abort("Invalid family")
    )
  }
  family <- switch(
    phenotype_class,
    "binary" = "binomial",
    "survival" = "cox",
    "continuous" = "gaussian"
  )

  # -- input validation -----------------------------------------------------
  chk::chk_is(matched_bulk, c("matrix", "data.frame"))
  chk::chk_is(sc_data, "Seurat")
  chk::chk_character(label_type)
  chk::chk_range(cutoff)
  chk::chk_list(reliability_test)
  chk::chk_list(cell_evaluation)

  # -- process dots ---------------------------------------------------------
  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption("verbose")
  seed <- dots$seed %||% SigBridgeRUtils::getFuncOption("seed")
  assay <- dots$assay %||% "RNA"
  load_cache <- dots$load_cache %||%
    dirname(dots$path2load_scissor_cache)
  save_cache <- dots$save_cache %||%
    dirname(dots$path2save_scissor_inputs)

  # -- fill defaults for reliability_test & cell_evaluation -----------------
  reliability_test <- utils::modifyList(
    list(run = FALSE, n = 10L, nfold = 10L),
    reliability_test
  )
  cell_evaluation <- utils::modifyList(
    list(
      run = FALSE,
      benchmark_data = "path_to_file.RData",
      FDR = 0.05,
      bootstrap_n = 100L
    ),
    cell_evaluation
  )

  # -- resolve label_type -------------------------------------------
  label_type_scissor <- if (family %chin% c("binomial", "cox")) {
    c(
      glue::glue("{label_type}_Negative"),
      glue::glue("{label_type}_Positive")
    )
  } else if (family == "gaussian") {
    n <- length(table(phenotype))
    glue::glue("{label_type}_{seq_len(n)}")
  }

  cache_config <- ScreenMethodConfig(
    method_name = "Scissor",
    param = get_env_vars(exclude = c("matched_bulk", "sc_data", "phenotype")),
    phenotype_class = phenotype_class,
    label_type = label_type
  )

  get_env_vars(exclude = c("matched_bulk", "sc_data", "phenotype")) # contains `cache_config`
}

#' @title Perform Scissor Screening Analysis
#' @description
#' Identifies phenotype-associated cell subpopulations in single-cell data using
#' regularized regression on matched bulk expression profiles. Scissor integrates
#' bulk and single-cell RNA-seq data to identify cells that are significantly
#' associated with phenotypic outcomes.
#'
#'
#' @inheritParams Screen
#' @param alpha Parameter used to balance the effect of the l1 norm and the network-based penalties. It can be a number or a searching vector. If alpha = NULL, a default searching vector is used. The range of alpha is between 0 and 1. A larger alpha lays more emphasis on the l1 norm.
#' @param cutoff  (default: `0.2`). When `alpha=NULL`, the cutoff is used to determine the optimal alpha.
#'        Higher values increase specificity.
#' @param family Model family for outcome type:
#'        - "gaussian": Continuous outcomes
#'        - "binomial": Binary outcomes (default)
#'        - "cox": Survival outcomes
#' @param reliability_test List controlling reliability testing:
#' - run: Whether to perform reliability test (default: FALSE)
#' - n: Permutation times (default: `10L`)
#' - nfold: Cross-validation folds (default: `10L`)
#' @param cell_evaluation List controlling cell evaluation:
#' - run: Whether to perform cell evaluation (default: FALSE)
#' - benchmark_data: Path to benchmark data (RData file)
#' - FDR_cutoff: FDR threshold for evaluation (default: `0.05`)
#' - bootstrap_n: Bootstrap iterations (default: `100L`)
#' @param ... Additional arguments. Currently supports:
#'    - `verbose`: Logical indicating whether to print progress messages. Defaults to `TRUE`.
#'    - `seed`: For reproducibility, default is `123L`
#'    - `assay`: Assay to use for single-cell data. Defaults to `"RNA"`
#'    - `load_cache`: Cache directory path for loading precomputed Scissor inputs.
#'      Supports root-level, cache-level, or parent-level paths. See [CacheSetHere()].
#'    - `save_cache`: Cache directory path for saving Scissor inputs. Supports
#'      root-level or parent-level paths. See [CacheSetHere()].
#'    - `path2load_scissor_cache` / `path2save_scissor_inputs`: Deprecated names,
#'      kept for backward compatibility.
#'
#' @return A list containing:
#' \describe{
#'   \item{scRNA_data}{A Seurat object with screened cells containing metadata:
#'     \describe{
#'       \item{scissor}{"Positive"/"Negative"/"Neutral" classification}
#'       \item{label_type}{Outcome label used}
#'     }
#'   }
#'   \item{scissor_result}{Raw Scissor results}
#'   \item{reliability_result}{If reliability_test=TRUE, contains:
#'     \describe{
#'       \item{statistic}{A value between 0 and 1}
#'       \item{p}{p-value of the test statistic}
#'       \item{AUC_test_real}{10 values of AUC for real data}
#'       \item{AUC_test_back}{A list of AUC for background data}
#'     }
#'   }
#'   \item{cell_evaluation}{If cell_evaluation=TRUE, contains:
#'     \describe{
#'       \item{evaluation_res}{A data.frame with some supporting information for each Scissor selected cell}
#'     }
#'   }
#' }
#'
#' @references
#' Sun D, Guan X, Moran AE, Wu LY, Qian DZ, Schedin P, et al. Identifying phenotype-associated subpopulations by integrating bulk and single-cell sequencing data. Nat Biotechnol. 2022 Apr;40(4):527–38.
#'
#' @section LICENSE:
#' Licensed under the GNU General Public License version 3 (GPL-3.0).
#' A copy of the license is available at <https://www.gnu.org/licenses/gpl-3.0.en.html>.
#'
#' @examples
#' \dontrun{
#' # Binary outcome example
#' res <- DoScissor(
#'   matched_bulk = bulk_matrix,
#'   sc_data = seurat_obj,
#'   phenotype = a_named_vector,
#'   family = "binomial"
#' )
#' }
#'
#' @export
#' @family screen_method
#' @family scissor
#'
DoScissor <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = "scissor",
  phenotype_class = c("binary", "continuous", "survival"),
  alpha = c(0.05, NULL),
  cutoff = 0.2,
  family = lifecycle::deprecated(),
  reliability_test = list(
    run = FALSE,
    n = 10L,
    nfold = 10L
  ),
  cell_evaluation = list(
    run = FALSE,
    benchmark_data = "path_to_file.RData",
    FDR_cutoff = 0.05,
    bootstrap_n = 100L
  ),
  path2load_scissor_cache = lifecycle::deprecated(),
  path2save_scissor_inputs = lifecycle::deprecated(),
  ...
) {
  # -- validate & prepare all parameters -----------------------------------
  p <- exec(ValidateScissorParams, !!!fn_fmls())

  # -- resolve cache paths using the caching system ------------------------
  load_file <- if (!is.null(p$load_cache)) {
    load_cache_dir <- CacheSetHere(
      path = p$load_cache,
      cache_config = p$cache_config,
      mode = "load",
      timestamp = p$dots$timestamp,
    )

    CheckCache(
      path = load_cache_dir,
      cache_config = p$cache_config
    )

    file.path(load_cache_dir, "Scissor_inputs.RData")
  } else {
    NULL
  }

  if (!is.null(p$save_cache)) {
    save_cache_dir <- CacheSetHere(
      path = p$save_cache,
      cache_config = p$cache_config,
      timestamp = p$dots$timestamp,
      mode = "save"
    )

    save_file <- file.path(save_cache_dir, "Scissor_inputs.RData")
  } else {
    save_file <- NULL
  }

  infos1 <- Scissor::Scissor.v5.optimized(
    bulk_dataset = matched_bulk,
    sc_dataset = sc_data,
    phenotype = phenotype,
    tag = p$label_type,
    alpha = alpha,
    cutoff = cutoff,
    family = p$family,
    Save_file = save_file,
    Load_file = load_file,
    verbose = p$verbose,
    seed = p$seed,
    assay = p$assay
  )

  # meta.data to add
  sc_meta <- data.frame(
    scissor = rep("Neutral", ncol(sc_data)),
    row.names = colnames(sc_data)
  )
  sc_meta$scissor[rownames(sc_meta) %chin% infos1$Scissor_pos] <- "Positive"
  sc_meta$scissor[rownames(sc_meta) %chin% infos1$Scissor_neg] <- "Negative"

  # * reliability test
  reliability_result <- if (p$reliability_test$run) {
    DoScissorRelTest(
      scissor_res = infos1,
      alpha = infos1$para$alpha,
      family = p$family,
      cell_num = length(infos1$Scissor_pos) +
        length(infos1$Scissor_neg),
      n = p$reliability_test$n,
      nfold = p$reliability_test$nfold,
      verbose = p$verbose
    )
  } else {
    NULL
  }

  # * cell_evaluation
  evaluate_res <- if (p$cell_evaluation$run) {
    DoScissorCellEval(
      benchmark_data_path = p$cell_evaluation$benchmark_data,
      scissor_res = infos1,
      FDR_cutoff = p$cell_evaluation$FDR_cutoff,
      bootstrap_n = p$cell_evaluation$bootstrap_n,
      verbose = p$verbose
    )
  } else {
    NULL
  }

  if (!is.null(p$save_cache)) {
    WriteCacheMeta(
      file = file.path(save_cache_dir, "cache_config.json"),
      cache_config = p$cache_config,
      additional_description = p$dots$additional_description,
      verbose = p$verbose,
    )
  }

  sc_data <- SeuratObject::AddMetaData(object = sc_data, metadata = sc_meta)
  sc_data <- AddMisc(
    seurat_obj = sc_data,
    Scissor = props(p$cache_config),
    cover = FALSE
  )

  list(
    scRNA_data = sc_data,
    scissor_result = infos1, # parameters included
    reliability_result = reliability_result,
    cell_evaluation = evaluate_res
  )
}

#' @title Perform Scissor Reliability Test
#' @description
#' Performs reliability testing for Scissor results to assess the stability
#' and robustness of identified phenotype-associated cells.
#'
#' @param scissor_res Scissor results from Scissor
#' @param alpha Alpha parameter used in Scissor, a scalar value between 0 and 1
#' @param family Model family used in Scissor, one of "gaussian", "binomial", "cox"
#' @param cell_num Number of cells identified by Scissor
#' @param n Permutation times (default: `10L`)
#' @param nfold Cross-validation folds (default: `10L`)
#' @param verbose Whether to show progress messages (default: `TRUE`)
#' @param ... No used parameters
#'
#' @return Reliability test results including statistics and p-values
#'
#' @keywords internal
#' @family scissor
DoScissorRelTest <- function(
  scissor_res,
  alpha,
  family,
  cell_num = length(scissor_res$Scissor_pos) +
    length(scissor_res$Scissor_neg),
  n = 10L,
  nfold = 10L,
  verbose = getFuncOption('verbose'),
  ...
) {
  skip_flag <- purrr::map_lgl(
    c(n, nfold),
    function(x) {
      if (!is.numeric(x) || is.na(x) || x != floor(x)) {
        cli::cli_warn(c(
          "x" = "`n` and `nfold` must be scalar integer. Skipping reliability test.",
          ">" = "Current `n`: {class(n)} {.val {n}}, `nfold`: {class(nfold)} {.val {nfold}}"
        ))
        return(TRUE)
      }
      FALSE
    }
  )
  skip_flag <- c(skip_flag, !chk::vld_range(alpha))
  if (any(skip_flag)) {
    cli::cli_warn("Skipping reliability test because of invalid parameters")
    return(NULL)
  }

  # indicate that Y has two levels, both Pos and Neg cells exist
  if (length(table(scissor_res$Y)) < 2) {
    cli::cli_warn(c(
      "x" = "Only one level detected in Scissor result. Skipping reliability test."
    ))
    return(NULL)
  }

  # * start
  if (verbose) {
    ts_cli$cli_alert_info(
      cli::col_green("Start reliability test")
    )
  }

  if (IsSkewedDynamic(scissor_res$Y)) {
    cli::cli_warn("Skewed distribution detected. Result may be unstable")
  }

  rel_res <- Scissor::reliability.test(
    X = scissor_res$X,
    Y = scissor_res$Y,
    network = scissor_res$network,
    alpha = alpha,
    family = family,
    cell_num = cell_num,
    n = n,
    nfold = nfold,
    verbose = verbose
  )

  if (verbose) {
    ts_cli$cli_alert_info(
      cli::col_green("reliability test: Done")
    )
  }

  rel_res
}

#' @title Perform Scissor Cell Evaluation
#' @description
#' Evaluates the significance of Scissor-identified cells using benchmark data
#' and bootstrap resampling.
#'
#' @param benchmark_data_path Path to benchmark data (RData file)
#' @param scissor_res Scissor results from Scissor
#' @param FDR_cutoff FDR threshold for significance (default: `0.05`)
#' @param bootstrap_n Number of bootstrap iterations (default: `100L`)
#' @param verbose Whether to show progress messages (default: `TRUE`)
#' @param ... No usage
#'
#' @return Cell evaluation results with significance assessments
#'
#' @keywords internal
#' @family scissor
DoScissorCellEval <- function(
  benchmark_data_path = 'path_to_file.RData',
  scissor_res,
  FDR_cutoff = 0.05,
  bootstrap_n = 100L,
  verbose = getFuncOption('verbose'),
  ...
) {
  if (!file.exists(benchmark_data_path)) {
    cli::cli_warn(c(
      'x' = '`benchmark_data` does not exist. Skipping cell evaluation'
    ))
    return(NULL)
  }
  if (FDR_cutoff <= 0 || FDR_cutoff >= 1) {
    cli::cli_warn(c(
      'x' = '`FDR_cutoff` must be between 0 and 1. Skipping cell evaluation',
      '>' = 'Current `FDR_cutoff`: {class(FDR_cutoff)} {.val {FDR_cutoff}}'
    ))
    return(NULL)
  }
  if (
    !is.numeric(bootstrap_n) ||
      is.na(bootstrap_n) ||
      bootstrap_n != floor(bootstrap_n)
  ) {
    cli::cli_warn(c(
      'x' = '`bootstrap_n` must be scalar integer. Skipping cell evaluation',
      '>' = 'Current `bootstrap_n`: {class(bootstrap_n)} {.val {bootstrap_n}}'
    ))
    return(NULL)
  }

  # * start
  if (verbose) {
    ts_cli$cli_alert_info(
      cli::col_green("Start cell evalutaion")
    )
  }

  evaluate_res <- Scissor::evaluate.cell(
    Load_file = benchmark_data_path,
    Scissor_result = scissor_res,
    FDR_cutoff = FDR_cutoff,
    bootstrap_n = bootstrap_n
  )

  if (verbose) {
    ts_cli$cli_alert_info(
      cli::col_green("Cell evalutaion finished")
    )
  }

  evaluate_res
}
