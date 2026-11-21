#' @title Benchmark Screening Methods: Classification Accuracy
#'
#' @description
#' Computes binary classification metrics for every screening method present in
#' the input, by comparing each method's per-cell labels against ground-truth
#' labels. `"Positive"` is treated as the positive class and all other labels
#' as negative, so the accuracy of different screening methods can be compared
#' on the same cells.
#'
#' @param obj The object to evaluate, one of:
#'   * a `Seurat` object whose `meta.data` contains screening method results
#'     (columns named after screening methods, e.g. `"Scissor"`);
#'   * a [ScreenMethodResult] object, whose `scRNA_data` slot is evaluated;
#'   * a [BenchAccConfig] object, whose `screen_method` entries are all run
#'     through [Screen()] first and the accumulated labels are evaluated in a
#'     single pass.
#' @param bench_config `character`. Ground-truth labels used as the reference:
#'   either a length-one column name in `obj@meta.data`, or a per-cell label
#'   vector of length `ncol(obj)`. Must contain a `"Positive"` label. For a
#'   [BenchAccConfig] input the labels are taken from the same-named slot.
#' @param verbose Logical. Whether to report progress information.
#'   Default: `TRUE`.
#' @param nthreads Integer. Number of threads used by the C++ metric
#'   computation. Default: `2L`.
#' @param ... Additional arguments passed between the S3 methods. Currently
#'   unused apart from forwarding
#'
#' @return A named `list` with one element per screening method found in
#'   `meta.data`. Each element is a 1-row numeric matrix with columns `TPR`,
#'   `FPR`, `TNR`, `FNR`, `F1`, `Accuracy` and so on.
#'
#' @details
#' ## Dispatch
#' `BenchAcc()` dispatches on the class of `obj`: a [BenchAccConfig] is screened
#' first, a [ScreenMethodResult] is unwrapped to its `scRNA_data` slot, and
#' `BenchAcc.Seurat()` performs the computation. Unsupported classes abort with
#' an informative error.
#'
#' ## Method detection
#' Columns of `obj@meta.data` whose names match a supported screening method
#' (case-insensitive) are taken as the predicted labels; other columns are
#' ignored. The function aborts when no such column exists.
#'
#' ## Ground truth
#' A length-one `bench_config` is looked up as a column of `obj@meta.data`;
#' otherwise it is used directly as a per-cell label vector. The labels must
#' include `"Positive"`, and the vector form must have one element per cell.
#'
#' ## Metrics
#' With `P`, `N`, `TP`, `FP`, `TN` and `FN` denoting the numbers of positive
#' cells, negative cells, true/false positives and true/false negatives
#'
#' @examplesIf rlang::is_installed("Seurat")
#' counts <- matrix(
#'   rep(1:10, 6L),
#'   nrow = 10L,
#'   dimnames = list(paste0("Gene", 1:10), paste0("Cell", 1:6))
#' )
#' obj <- Seurat::CreateSeuratObject(counts = counts)
#' obj$truth <- c("Positive", "Negative", "Positive", "Other", "Positive", "Negative")
#' obj$Scissor <- c("Positive", "Negative", "Positive", "Other", "Negative", "Negative")
#' BenchAcc(obj, bench_config = "truth", verbose = FALSE)
#'
#' @seealso [BenchAccConfig()], [ScreenMethodResult], [Screen()],
#'   [BenchGSEA()], [BenchMemory()], [BenchGeneScore()]
#'
#' @family Benchmark
#'
#' @name BenchAcc
#' @export
BenchAcc <- function(
  obj,
  ...
) {
  UseMethod("BenchAcc")
}

#' @rdname BenchAcc
#' @export
BenchAcc.Seurat <- function(
  obj,
  bench_config = character(),
  ...,
  verbose = getFuncOption("verbose"),
  nthreads = getFuncOption("nthreads")
) {
  meta_data <- obj[[]]
  chk::chk_character(bench_config)
  if (length(bench_config) == 0L) {
    Abort(
      "{.val bench_config} is not specified"
    )
  } else if (length(bench_config) == 1L) {
    if (!bench_config %chin% colnames(meta_data)) {
      Abort(
        "{.val {bench_config}} is not a valid column name in obj@meta.data"
      )
    }
    bench_label <- meta_data[[bench_config]]
    bench_label_names <- names(table(bench_label))
    if (!"Positive" %chin% bench_label_names) {
      Abort(
        "{.val {bench_config}} in obj@meta.data doesn't have a `Positive` label"
      )
    }
  } else {
    bench_label <- bench_config
    bench_label_names <- names(table(bench_label))

    if (!"Positive" %chin% bench_label_names) {
      Abort(
        "{.val bench_config} doesn't have a `Positive` label"
      )
    }
    if (length(bench_label) != ncol(obj)) {
      n_cells <- ncol(obj)
      n_bench_label <- length(bench_label)
      Abort(
        "{.val bench_config} doesn't have the same length as cells",
        "Number of cells: {.val {n_cells}}, length of bench_config: {.val {n_bench_label}}"
      )
    }
  }
  available_method <- names(ScreenStrategy)
  existing_method <- tolower(colnames(meta_data)) %chin% tolower(available_method)
  if (!any(existing_method)) {
    Abort(
      "Seurat object has no screen method results",
      "Can't find any column names like {.val {available_method}}"
    )
  }
  metrics <- meta_data[, existing_method]

  if (verbose) {
    col_names <- colnames(metrics)
    cli::cli_alert_info("Calculating metrics for {.val {col_names}}")
  }

  # cpp
  calc_binary_metrics(
    metrics = metrics,
    bench_label = bench_label,
    n_threads = nthreads
  )
  # a list of matrices of metrics
}

#' @rdname BenchAcc
#' @export
`BenchAcc.SigBridgeR::BenchAccConfig` <- function(
  obj,
  ...,
  verbose = getFuncOption("verbose")
) {
  sc <- obj@sc_data

  for (method in obj@screen_method) {
    if (verbose) {
      cli::cli_progress_step(
        msg = "Fetching {.val {method}} results",
        msg_done = "Fetched {.val {method}} results",
        msg_failed = "Failed to fetch {.val {method}} results"
      )
    }
    screen_res <- Screen(
      matched_bulk = obj@matched_bulk,
      sc_data = sc,
      phenotype = obj@phenotype,
      label_type = obj@label_type,
      phenotype_class = obj@phenotype_class,
      screen_method = method,
      !!!obj@method_param[[paste0(method, "_param")]]
    )
    sc <- screen_res@scRNA_data
    gc(verbose = FALSE)
  }
  if (verbose) {
    cli::cli_progress_done()
  }

  BenchAcc.Seurat(
    obj = sc,
    bench_config = obj@bench_config, # yes it's a chr vec
    ...,
    verbose = verbose
  )
}

#' @rdname BenchAcc
#' @export
`BenchAcc.SigBridgeR::ScreenMethodResult` <- function(obj, ...) {
  BenchAcc.Seurat(obj = obj@scRNA_data, ...)
}

#' @rdname BenchAcc
#' @export
BenchAcc.default <- function(obj, ...) {
  cls_obj <- class(obj)
  avaialble_methods <- utils::methods(BenchAcc, all.names = TRUE) |>
    as.character() |>
    gsub(pattern = ".*\\.", replacement = "", x = _)
  Abort(
    "BenchAcc is not implemented for class: {.cls {cls_obj}}",
    "Supported classes: {.cls {avaialble_methods}}"
  )
}
