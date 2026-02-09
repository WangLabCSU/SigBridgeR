# ? ---- 2. DO SCISSOR ----

#' @title Perform Scissor Screening Analysis
#' @description
#' Identifies phenotype-associated cell subpopulations in single-cell data using
#' regularized regression on matched bulk expression profiles. Scissor integrates
#' bulk and single-cell RNA-seq data to identify cells that are significantly
#' associated with phenotypic outcomes.
#'
#' @usage
#' DoScissor(
#'    path2load_scissor_cache = NULL,
#'    path2save_scissor_inputs = "Scissor_inputs.RData",
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
#' @param path2load_scissor_cache Path to precomputed Scissor inputs (RData file).
#'        If provided, skips recomputation (default: NULL).
#' @param path2save_scissor_inputs Path to save intermediate files (default: "Scissor_inputs.RData").
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
#'    - `assay`: Assay to use for single-cell data. Defaults to `"RNA"
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
  path2load_scissor_cache = NULL,
  path2save_scissor_inputs = "Scissor_inputs.RData",
  matched_bulk,
  sc_data,
  phenotype,
  label_type = "scissor",
  alpha = c(0.05, NULL),
  cutoff = 0.2,
  family = c("gaussian", "binomial", "cox"),
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
  ...
) {
  # * Input validation
  chk::chk_is(matched_bulk, c("matrix", "data.frame"))
  chk::chk_is(sc_data, "Seurat")
  chk::chk_character(label_type)
  chk::chk_range(cutoff)
  family <- SigBridgeRUtils::MatchArg(
    family,
    c("gaussian", "binomial", "cox"),
    NULL
  )
  chk::chk_list(reliability_test)
  chk::chk_list(cell_evaluation)

  # * get defaults from dots
  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption("verbose")
  seed <- dots$seed %||% SigBridgeRUtils::getFuncOption("seed")
  assay <- dots$assay %||% "RNA"

  # * default setting for `reliability_test` & `cell_evaluation`
  default_reliability_test <- list(
    run = FALSE,
    n = 10L,
    nfold = 10L
  )
  default_cell_evalutaion <- list(
    run = FALSE,
    benchmark_data = "path_to_file.RData",
    FDR = 0.05,
    bootstrap_n = 100L
  )
  reliability_test <- utils::modifyList(
    default_reliability_test,
    reliability_test
  )
  cell_evaluation <- utils::modifyList(
    default_cell_evalutaion,
    cell_evaluation
  )

  if (family %chin% c("binomial", "cox")) {
    label_type_scissor <- c(
      glue::glue("{label_type}_Negative"),
      glue::glue("{label_type}_Positive")
    )
  } else if (family == "gaussian") {
    n <- length(table(phenotype))
    label_type_scissor <- glue::glue("{label_type}_{seq_len(n)}")
  }

  if (!is.null(path2save_scissor_inputs)) {
    path <- dirname(path2save_scissor_inputs)
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
    }
  }

  infos1 <- Scissor::Scissor.v5.optimized(
    bulk_dataset = matched_bulk,
    sc_dataset = sc_data,
    phenotype = phenotype,
    tag = label_type_scissor,
    alpha = alpha,
    cutoff = cutoff,
    family = family,
    Save_file = path2save_scissor_inputs,
    Load_file = path2load_scissor_cache,
    verbose = verbose,
    seed = seed,
    assay = assay
  )

  # meta.data to add
  sc_meta <- data.frame(
    scissor = rep("Neutral", ncol(sc_data)),
    row.names = colnames(sc_data)
  )
  sc_meta$scissor[rownames(sc_meta) %chin% infos1$Scissor_pos] <- "Positive"
  sc_meta$scissor[rownames(sc_meta) %chin% infos1$Scissor_neg] <- "Negative"
  sc_data <- Seurat::AddMetaData(sc_data, metadata = sc_meta) %>%
    SigBridgeRUtils::AddMisc(
      scissor_type = label_type,
      scissor_para = infos1$para,
      cover = FALSE
    )

  # * reliability test
  reliability_result <- if (reliability_test$run) {
    DoScissorRelTest(
      scissor_res = infos1,
      alpha = alpha,
      family = family,
      cell_num = length(infos1$Scissor_pos) +
        length(infos1$Scissor_neg),
      n = reliability_test$n,
      nfold = reliability_test$nfold,
      verbose = verbose
    )
  } else {
    NULL
  }

  # * cell_evaluation
  evaluate_res <- if (cell_evaluation$run) {
    DoScissorCellEval(
      benchmark_data_path = cell_evaluation$benchmark_data,
      scissor_res = infos1,
      FDR_cutoff = cell_evaluation$FDR_cutoff,
      bootstrap_n = cell_evaluation$bootstrap_n,
      verbose = verbose
    )
  } else {
    NULL
  }

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
  if (any(skip_flag)) {
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
