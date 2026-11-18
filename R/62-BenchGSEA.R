#' @title Benchmark Screening Methods: Gene Set Enrichment Analysis
#'
#' @description
#' Runs gene set enrichment analysis (GSEA) on the marker genes of every
#' screening method present in the input, so that the functional signals
#' recovered by different methods can be compared side by side.
#'
#' @param obj A `Seurat` object whose `meta.data` contains screening method
#'   results (columns named after screening methods, e.g. `"Scissor"`), a
#'   `ScreenMethodResult` object, or a [BenchGSEAConfig] object.
#' @param group A vector of length equal to the number of cells, with each
#'   element being the group label of the corresponding cell.
#' @param pathway A named list of gene sets, each element a character vector of
#'   gene symbols. Passed to [fgsea::fgsea()].
#' @param p_value_threshold Numeric or `NULL`. When not `NULL`, markers are
#'   filtered by adjusted p-value (`p_val_adj < 0.05`) before ranking.
#'   Default: `NULL`.
#' @param verbose Logical. Whether to report progress.
#'   Default: `TRUE`.
#' @param ... Additional arguments, automatically split by formal argument
#'   names between [Seurat::FindMarkers()] and [fgsea::fgsea()] (including
#'   [fgsea::fgseaMultilevel()] / [fgsea::fgseaSimple()]).
#'   `ident.1`, `ident.2` and `group.by` are reserved by `BenchGSEA` and must
#'   not be supplied.
#'
#' @return A `data.table` of [fgsea::fgsea()] results stacked across screening
#'   methods, with an extra `screenMethod` column identifying the method and a
#'   `leadingEdgeStats` list column holding the ranking statistics of the
#'   leading-edge genes. Methods without markers are dropped with a warning.
#'
#' @details
#' For each screening method found in `meta.data`:
#'
#' 1. Labels are binarised to `"Positive"` vs `"Other"`;
#' 2. [Seurat::FindMarkers()] is run with `ident.1 = "Positive"`,
#'    `ident.2 = "Other"` and `group.by` set to the method column;
#' 3. Markers are ranked by `sign(avg_log2FC) * -log10(p_val_adj)`;
#' 4. [fgsea::fgsea()] is run on the ranked list, and the ranking statistics of
#'    the leading-edge genes are attached as `leadingEdgeStats`.
#'
#' For a [BenchGSEAConfig] input, every configured screening method is first
#' run via [Screen()] on the same data before the GSEA comparison.
#'
#' @examples
#' \dontrun{
#' # sc is a Seurat object already screened by e.g. Scissor and scPAS
#' pathways <- list(
#'   pathway_a = c("Gene1", "Gene2", "Gene3"),
#'   pathway_b = c("Gene4", "Gene5", "Gene6")
#' )
#' res <- BenchGSEA(sc, pathway = pathways)
#' }
#'
#' @seealso [BenchGSEAConfig()], [Screen()], [fgsea::fgsea()],
#'   [Seurat::FindMarkers()], [BenchAcc()], [BenchMemory()]
#'
#' @family Benchmark
#'
#' @export
#'
#' @name BenchGSEA
BenchGSEA <- function(obj, ...) {
  UseMethod("BenchGSEA")
}

#' @rdname BenchGSEA
#' @export
BenchGSEA.Seurat <- function(
  obj,
  group,
  pathway,
  p_value_threshold = NULL,
  ...,
  verbose = getFuncOption("verbose")
) {
  rlang::check_installed(c("dplyr", "fgsea"))

  dots <- list(...)
  dots_names <- names(dots)

  fgsea_fmls_names <- union(rlang::fn_fmls_names(fgsea::fgseaMultilevel), rlang::fn_fmls_names(fgsea::fgseaSimple))
  fgsea_params <- dots[dots_names %chin% fgsea_fmls_names]
  find_markers_params <- dots[!dots_names %chin% fgsea_fmls_names]

  if (any(c("ident.1", "ident.2", "group.by") %chin% dots_names)) {
    Abort(
      "Can't use {.arg ident.1} or {.arg ident.2} or {.arg group.by}, they are reserved for BenchGSEA",
      "Please remove {.arg ident.1} or {.arg ident.2} or {.arg group.by}"
    )
  }

  # ------------------------------------------------------------------------------------------------------------------------------------

  meta_data <- obj[[]]
  existing_method <- tolower(colnames(meta_data)) %chin% tolower(AVAILABLE_SCREEN_METHOD)
  if (!any(existing_method)) {
    Abort(
      "Seurat object has no screen method results",
      "Can't find any column names like {.val {AVAILABLE_SCREEN_METHOD}}"
    )
  }
  # * binarize
  existing_method_names <- colnames(meta_data)[existing_method]
  metrics <- meta_data |>
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::all_of(existing_method_names),
        .fns = function(label) {
          dplyr::case_when(
            label != "Positive" ~ "Other",
            label == "Positive" ~ "Positive"
          )
        }
      )
    )
  obj[[]] <- metrics # tmp modify

  v <- NULL # for progress bar
  markers <- lapply(X = existing_method_names, FUN = function(method) {
    if (verbose) {
      cli::cli_progress_step(
        "Finding markers for {.val {method}} groups",
        "Found {.val {length(v)}} markers for {.val {method}} groups"
      )
    }
    markers <- do.call(
      what = Seurat::FindMarkers,
      args = rlang::list2(
        object = obj,
        ident.1 = "Positive",
        ident.2 = "Other",
        group.by = method,
        !!!find_markers_params
      )
    )
    if (!is.null(p_value_threshold)) {
      markers <- markers[markers$p_val_adj < 0.05, ]
    }

    v <- setNames(
      object = sign(markers$avg_log2FC) * -log10(markers$p_val_adj),
      nm = rownames(markers)
    )

    if (verbose) {
      cli::cli_progress_update()
    }

    list(filtered = v, all = markers)
  })
  names(markers) <- existing_method_names

  if (verbose) {
    cli::cli_progress_step("Calculating GSEA", "Done")
  }

  leadingEdgeStats <- leadingEdge <- NULL # ease checking NOTE
  fgsea_res <- purrr::map(X = markers, FUN = function(markers, name) {
    if (length(markers$filtered) == 0L) {
      cli::cli_warn("No markers found for {name}")
      return(NULL)
    }

    dt <- do.call(
      what = fgsea::fgsea,
      args = list(pathways = pathway, stats = markers$filtered, !!!fgsea_params)
    )
    dt[,
      leadingEdgeStats := lapply(
        X = leadingEdge,
        FUN = function(genes) markers$filtered[genes]
      )
    ]
    dt
  })

  if (verbose) {
    cli::cli_progress_done()
  }

  data.table::rbindlist(fgsea_res, idcol = "screenMethod")
}


#' @rdname BenchGSEA
#' @export
`BenchGSEA.SigBridgeR::ScreenMethodResult` <- function(obj, ...) {
  BenchGSEA.Seurat(obj = obj, ...)
}

#' @rdname BenchGSEA
#' @export
`BenchGSEA.SigBridgeR::BenchGSEAConfig` <- function(
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

  BenchGSEA.Seurat(
    obj = sc,
    ...,
    verbose = verbose
  )
}

#' @rdname BenchGSEA
#' @export
BenchGSEA.default <- function(obj, ...) {
  cls_obj <- class(obj)
  avaialble_methods <- utils::methods(BenchGSEA, all.names = TRUE) |>
    as.character() |>
    gsub(pattern = ".*\\.", replacement = "", x = _)
  Abort(
    "BenchGSEA is not implemented for class: {.cls {cls_obj}}",
    "Available methods: {.cls {avaialble_methods}}"
  )
}
