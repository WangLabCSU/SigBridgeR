#' @return A list containing metrics matrix, with column `TPR FPR TNR FNR F1 Accuracy`
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
  rlang::check_string(bench_config)
  if (length(bench_config) == 0L) {
    Abort(
      "{.val bench_config} is not specified"
    )
  } else if (length(bench_config) == 1L) {
    if (!bench_config %chin% colnames(meta_data)) {
      Abort(
        "{.val {bench_config}} is not a valid column name in sc_data@meta.data"
      )
    }
    bench_label <- meta_data[[bench_config]]
    bench_label_names <- names(table(bench_label))
    if (!"Positive" %chin% bench_label_names) {
      Abort(
        "{.val {bench_config}} in sc_data@meta.data doesn't have a `Positive` label"
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
    if (length(bench_label) != ncol(sc_data)) {
      n_cells <- ncol(sc_data)
      n_bench_label <- length(bench_label)
      Abort(
        "{.val bench_config} doesn't have the same length as cells",
        "Number of cells: {.val {n_cells}}, length of bench_config: {.val {n_bench_label}}"
      )
    }
  }
  existing_method <- tolower(colnames(meta_data)) %chin% tolower(AVAILABLE_SCREEN_METHOD)
  if (!any(existing_method)) {
    Abort(
      "Seurat object has no screen method results",
      "Can't find any column names like {.val {AVAILABLE_SCREEN_METHOD}}"
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
