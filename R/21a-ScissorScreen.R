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
#'   `label_type_scissor`, `phenotype_class_cache`, `dots`.
#'
#' @keywords internal
#' @family scissor
ValidateScissorParams <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type,
  phenotype_class,
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
      "DoScissor(family = c('gaussian', 'binomial', 'cox'))",
      "DoScissor(phenotype_class = c('binary', 'continuous', 'survival'))"
    )
    family <- SigBridgeRUtils::MatchArg(
      family,
      c("gaussian", "binomial", "cox"),
      NULL
    )
    phenotype_class <- family
  }

  # -- input validation -----------------------------------------------------
  chk::chk_is(matched_bulk, c("matrix", "data.frame"))
  chk::chk_is(sc_data, "Seurat")
  chk::chk_character(label_type)
  chk::chk_range(cutoff)
  chk::chk_list(reliability_test)
  chk::chk_list(cell_evaluation)

  phenotype_class <- SigBridgeRUtils::MatchArg(
    phenotype_class,
    c('binary', 'continuous', 'survival')
  )

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

  # -- resolve label_type_scissor -------------------------------------------
  if (phenotype_class %chin% c("binary", "survival")) {
    label_type_scissor <- c(
      glue::glue("{label_type}_Negative"),
      glue::glue("{label_type}_Positive")
    )
  } else {
    n <- length(table(phenotype))
    label_type_scissor <- glue::glue("{label_type}_{seq_len(n)}")
  }

  # -- map family -> phenotype_class ----------------------------------------

  list(
    family = family,
    verbose = verbose,
    seed = seed,
    assay = assay,
    load_cache = load_cache,
    save_cache = save_cache,
    reliability_test = reliability_test,
    cell_evaluation = cell_evaluation,
    label_type = label_type,
    label_type_scissor = label_type_scissor,
    phenotype_class = phenotype_class,
    dots = dots
  )
}

#' @title Perform Scissor Screening Analysis
#' @description
#' Identifies phenotype-associated cell subpopulations in single-cell data using
#' regularized regression on matched bulk expression profiles. Scissor integrates
#' bulk and single-cell RNA-seq data to identify cells that are significantly
#' associated with phenotypic outcomes.
#'
#' @usage
#' DoScissor(
#'    matched_bulk,
#'    sc_data,
#'    phenotype,
#'    label_type = "scissor",
#'    alpha = c(0.05, NULL),
#'    cutoff = 0.2,
#'    family = c("gaussian", "binomial", "cox"),
#'    reliability_test = list(
#'      run = FALSE, # whether to run reliability test
#'      n = 10L, # permutation times
#'      nfold = 10L # cross validation folds
#'    ),
#'    cell_evaluation = list(
#'      run = FALSE, # whether to run cell evaluation
#'      benchmark_data = "path_to_file.RData", # path to benchmark data
#'      FDR_cutoff = 0.05,
#'      bootstrap_n = 100L
#'    ),
#'    ...
#' )
#'
#' @param matched_bulk Normalized bulk expression matrix (features × samples).
#'        Column names must match `phenotype` identifiers.
#' @param sc_data Seurat object containing single-cell RNA-seq data.
#' @param phenotype Clinical outcome data. Can be:
#'        - Vector: named with sample IDs
#'        - Data frame: with row names matching bulk columns
#' @param label_type Character specifying phenotype label type (e.g., "SBS1", "time"), stored in `scRNA_data@misc`
#' @param alpha Parameter used to balance the effect of the l1 norm and the network-based penalties. It can be a number or a searching vector. If alpha = NULL, a default searching vector is used. The range of alpha is between 0 and 1. A larger alpha lays more emphasis on the l1 norm.
#' @param cutoff  (default: `0.2`). When `alpha=NULL`, the cutoff is used to determine the optimal alpha.
#'        Higher values increase specificity.
#' @param family `r lifecycle::badge('deprecated')` Model family for outcome type:
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
#'    - `r lifecycle::badge('deprecated')` `path2load_scissor_cache` / `path2save_scissor_inputs`: Deprecated names,
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
  phenotype_class = c("binary", "continuous", "survival"),
  path2load_scissor_cache = lifecycle::deprecated(),
  path2save_scissor_inputs = lifecycle::deprecated(),
  ...
) {
  # -- validate & prepare all parameters -----------------------------------
  p <- ValidateScissorParams(
    matched_bulk = matched_bulk,
    sc_data = sc_data,
    phenotype = phenotype,
    label_type = label_type,
    phenotype_class = phenotype_class,
    alpha = alpha,
    cutoff = cutoff,
    family = family,
    reliability_test = reliability_test,
    cell_evaluation = cell_evaluation,
    path2load_scissor_cache = path2load_scissor_cache,
    path2save_scissor_inputs = path2save_scissor_inputs,
    ...
  )

  if (!is.null(p$load_cache) || !is.null(p$save_cache)) {
    cache_config <- ScreenMethodConfig(
      method_name = "Scissor",
      phenotype_class = p$phenotype_class,
      label_type = p$label_type,
      param = fn_fmls()
    )
  }

  # -- resolve cache paths using the caching system ------------------------
  load_file <- if (!is.null(p$load_cache)) {
    load_cache_dir <- CacheSetHere(
      path = p$load_cache,
      cache_config = cache_config,
      timestamp = p$dots$timestamp,
      mode = "load"
    )

    CheckCache(
      path = load_cache_dir,
      cache_config = cache_config
    )

    file.path(load_cache_dir, "Scissor_inputs.RData")
  } else {
    NULL
  }

  save_file <- if (!is.null(p$save_cache)) {
    save_cache_dir <- CacheSetHere(
      path = p$load_cache,
      cache_config = cache_config,
      timestamp = p$dots$timestamp,
      mode = "save"
    )

    file.path(save_cache_dir, "Scissor_inputs.RData")
  } else {
    NULL
  }

  infos1 <- Scissor::Scissor.v5.optimized(
    bulk_dataset = matched_bulk,
    sc_dataset = sc_data,
    phenotype = phenotype,
    tag = p$label_type_scissor,
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
      cache_config = cache_config,
      verbose = p$verbose,
      additional_description = p$dots$additional_description
    )
  }

  sc_data <- SeuratObject::AddMetaData(object = sc_data, metadata = sc_meta)
  sc_data <- AddMisc(
    seurat_obj = sc_data,
    scissor_type = label_type,
    scissor_para = c(infos1$para, reliability_test = reliability_result),
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

  if (is_skewed_dynamic(scissor_res$Y)) {
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

#' @title Is 1 Much More Than 0
#' @description
#' A function to check whether Y is highly skewed. The expected ideal distribution is p1 : p0 = cutoff : (1 - cutoff).
#' If the deviation exceeds n_sd times the standard deviation, it is considered skewed,
#' i.e., the input phenotype raw data is skewed, which means the reliability test may be unreliable or error-prone.
#'
#' @keywords internal
is_skewed_dynamic <- function(
  x, # a vector
  target = 0, # background value
  expected_p = 0.8, # = cutoff
  n_sd = 4L # n times standard deviation
) {
  n <- length(x)
  if (n == 0) {
    return(NA)
  }

  p_hat <- mean(x == target)
  sd_expected <- sqrt(expected_p * (1 - expected_p) / n)

  # 偏离超过 n_sd 倍标准差则判定为偏态
  abs(p_hat - expected_p) > n_sd * sd_expected
}
