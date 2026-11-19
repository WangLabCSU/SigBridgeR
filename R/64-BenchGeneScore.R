# --------------------------------------------------------------------------------
# BenchGeneScore
# --------------------------------------------------------------------------------

#' @title Single-sample Gene Set Scoring for Single-Cell Data
#'
#' @description
#' Scores one or more gene sets in every cell of a screened single-cell dataset
#' with `ssGSEA`, `AUCell`, `UCell` or `AddModuleScore`. All backends return the
#' same gene-set-by-cell matrix layout, so the cell-level signal recovered by
#' different screening methods can be compared on a common ground.
#'
#' @param obj The object to score, one of:
#'   * a `Seurat` object, scored as is;
#'   * a [ScreenMethodResult] object, whose `scRNA_data` slot is scored;
#'   * a [BenchGeneScoreConfig] object, whose `screen_method` entries are all
#'     run through [Screen()] first and the accumulated labels are scored in a
#'     single pass.
#' @param pathway A non-empty, named list of character vectors, each element
#'   being the genes of one gene set. List names are used as the gene-set names
#'   (row names) of the result.
#' @param method Character. Scoring backend, one of `"ssGSEA"`, `"AUCell"`,
#'   `"UCell"` or `"AddModuleScore"`. Default: `"ssGSEA"`.
#' @param assay Character or `NULL`. Assay to score. Default: `NULL`, i.e. the
#'   default assay of `obj` as reported by
#'   [SeuratObject::DefaultAssay()].
#' @param layer Character. Expression layer of `assay` holding the values to
#'   score. Default: `"data"`.
#' @param verbose Logical. Whether to report progress information.
#'   Default: `TRUE`.
#' @param ... Additional named arguments forwarded to the selected backend.
#'   Arguments are matched against the formal argument names of the backend and
#'   split automatically over the two stages when a method needs two calls:
#'   * `"ssGSEA"`: [GSVA::ssgseaParam()] and [GSVA::gsva()];
#'   * `"AUCell"`: [AUCell::AUCell_buildRankings()] and
#'     [AUCell::AUCell_calcAUC()];
#'   * `"UCell"`: [UCell::ScoreSignatures_UCell()];
#'   * `"AddModuleScore"`: [Seurat::AddModuleScore()].
#'
#' @return A numeric matrix with gene sets in rows and cells in columns. Gene
#'   sets dropped by backend-specific filtering are omitted from the result.
#'
#' @details
#' ## Dispatch
#' `BenchGeneScore()` dispatches on the class of `obj`, so a
#' [BenchGeneScoreConfig] is benchmarked and a [ScreenMethodResult] is unwrapped
#' before the actual scoring happens. [BenchGeneScore.Seurat()] is itself a
#' generic dispatching on `method`: every backend lives in its own
#' `BenchGeneScoreSeurat.<method>` function and its package is loaded on demand
#' by [rlang::check_installed()], so neither `GSVA`, `AUCell`, `UCell` nor
#' `Seurat` is required to merely call the generic.
#'
#' ## Config input
#' For a [BenchGeneScoreConfig] input, the configured `screen_method` entries
#' are run one after another with [Screen()]. Every run adds its labels to the
#' `meta.data` of the working `Seurat` object and the labels are accumulated
#' across methods, after which the scores are computed once. Scores are
#' therefore comparable across screening methods by construction.
#'
#' ## Backend notes
#' * `"UCell"` is called with `name = ""` so that gene-set names stay untouched,
#'   and its cell-by-gene-set output is transposed to the common layout.
#' * `"AddModuleScore"` drops the genes absent from `assay`, drops the gene sets
#'   that become empty, and stores the scores under an internal
#'   `.BenchGeneScore` prefix so that the `meta.data` of `obj` is left untouched.
#'   It aborts when no gene set survives the filtering.
#' * `"AUCell"` builds the rankings with `plotStats = FALSE`, since the plot is
#'   only meaningful for a manual inspection.
#'
#' @examplesIf rlang::is_installed(c("Seurat", "GSVA"))
#' counts <- matrix(
#'   rep(seq_len(20), times = 5) + rep(seq_len(5), each = 20),
#'   nrow = 20,
#'   dimnames = list(paste0("Gene", seq_len(20)), paste0("Cell", seq_len(5)))
#' )
#' obj <- Seurat::CreateSeuratObject(counts = counts)
#' obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#' pathways <- list(
#'   First = paste0("Gene", 1:5),
#'   Second = paste0("Gene", 6:10)
#' )
#' scores <- BenchGeneScore(obj, pathways, method = "ssGSEA", verbose = FALSE)
#'
#' @seealso [BenchGeneScoreConfig()], [BenchGeneScoreConfigSet()],
#'   [Screen()], [BenchGSEA()], [BenchAcc()]
#'
#' @family Benchmark
#'
#' @name BenchGeneScore
#' @export
BenchGeneScore <- function(obj, ...) {
  UseMethod("BenchGeneScore")
}

#' @rdname BenchGeneScore
#' @export
`BenchGeneScore.SigBridgeR::BenchGeneScoreConfig` <- function(
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

  BenchGeneScore.Seurat(
    obj = sc,
    ...,
    verbose = verbose
  )
}

#' @rdname BenchGeneScore
#' @export
`BenchGeneScore.SigBridgeR::ScreenMethodResult` <- function(obj, ...) {
  BenchGeneScore.Seurat(obj = obj@scRNA_data, ...)
}

#' @rdname BenchGeneScore
#' @export
BenchGeneScore.Seurat <- function(
  obj,
  method = c("ssGSEA", "AUCell", "UCell", "AddModuleScore"),
  ...
) {
  method <- rlang::arg_match(method)
  switch(
    method,
    "ssGSEA" = {
      rlang::check_installed("GSVA")
      BenchGeneScoreSeurat.ssGSEA(
        obj = obj,
        ...
      )
    },
    "AUCell" = {
      rlang::check_installed("AUCell")
      BenchGeneScoreSeurat.AUCell(
        obj = obj,
        ...
      )
    },
    "UCell" = {
      rlang::check_installed("UCell")
      BenchGeneScoreSeurat.UCell(
        obj = obj,
        ...
      )
    },
    "AddModuleScore" = {
      BenchGeneScoreSeurat.AddModuleScore(
        obj = obj,
        ...
      )
    }
  )
}

# --------------------------------------------------------------------------------
# Backends
# --------------------------------------------------------------------------------

#' @title Gene Set Scoring Backends
#'
#' @description
#' Internal implementations behind [BenchGeneScore.Seurat()], one per supported
#' value of its `method` argument. They share the same interface and only differ
#' in the backend they wrap, hence the common documentation.
#'
#' `BenchGeneScoreSeurat.<method>()` is not meant to be called directly: it is
#' reached through the class that [BenchGeneScore.Seurat()] sets on `method`
#' before dispatching.
#'
#' @param obj A `Seurat` object.
#' @param ... Additional arguments forwarded to the wrapped backend.
#' @param verbose Logical. Whether to print progress information.
#'   Default: `getFuncOption("verbose")`.
#'
#' @return The score matrix of the wrapped backend, with gene sets in rows and
#'   cells in columns, except `"AUCell"`, which returns the
#'   [AUCell::AUCell_calcAUC()] object whose AUC matrix uses the same layout.
#'
#' @details
#' Expression values are always taken with [SeuratObject::LayerData()], so the
#' backend receives a plain gene-by-cell matrix. Each backend then does its own
#' housekeeping:
#' * `ssGSEA` splits `...` between [GSVA::ssgseaParam()] and [GSVA::gsva()];
#' * `AUCell` splits `...` between [AUCell::AUCell_buildRankings()] and
#'   [AUCell::AUCell_calcAUC()], and builds the rankings with
#'   `plotStats = FALSE`;
#' * `UCell` passes everything to [UCell::ScoreSignatures_UCell()] with
#'   `name = ""` so that gene-set names stay untouched, then transposes the
#'   cell-by-gene-set result;
#' * `AddModuleScore` keeps only the genes present in `assay`, drops the gene
#'   sets that become empty, and aborts when none of them is left.
#'
#' @seealso [BenchGeneScore()]
#'
#' @keywords internal
#' @name BenchGeneScoreSeuratMethod
BenchGeneScoreSeurat.ssGSEA <- function(
  obj,
  pathway,
  assay = NULL,
  layer = "data",
  ..., # pass to GSVA
  verbose = getFuncOption("verbose")
) {
  rlang::check_installed("GSVA", reason = "to compute ssGSEA scores")
  expr <- SeuratObject::LayerData(object = obj, assay = assay, layer = layer)

  dots <- list(...)

  ssgsea_param_names <- names(formals(GSVA::ssgseaParam))
  gsva_names <- names(formals(GSVA::gsva))

  ssgsea_param <- dots[names(dots) %chin% ssgsea_param_names]
  gsva_param <- dots[!names(dots) %chin% ssgsea_param_names] # receive dots

  param <- do.call(
    GSVA::ssgseaParam,
    modify_list(
      list(exprData = expr, geneSets = pathway),
      ssgsea_param
    )
  )
  scores <- do.call(
    GSVA::gsva,
    modify_list(
      list(param = param, verbose = verbose),
      gsva_param
    )
  )
  scores
}


#' @rdname BenchGeneScoreSeuratMethod
#' @keywords internal
BenchGeneScoreSeurat.AUCell <- function(
  obj,
  pathway,
  assay = NULL,
  layer = "data",
  ...,
  verbose = getFuncOption("verbose")
) {
  rlang::check_installed("AUCell", reason = "to compute AUCell scores")
  dots <- list(...)

  expr <- SeuratObject::LayerData(object = obj, assay = assay, layer = layer)

  build_rankings_names <- names(formals(AUCell::AUCell_buildRankings))
  calc_auc_names <- names(formals(AUCell::AUCell_calcAUC))

  calc_auc_param <- dots[!names(dots) %chin% build_rankings_names]
  build_rankings_param <- dots[names(dots) %chin% build_rankings_names]

  rankings <- do.call(
    AUCell::AUCell_buildRankings,
    modify_list(
      list(exprMat = expr, plotStats = FALSE, verbose = verbose),
      build_rankings_param
    )
  )
  auc <- do.call(
    AUCell::AUCell_calcAUC,
    modify_list(
      list(geneSets = pathway, rankings = rankings, verbose = verbose),
      calc_auc_param
    )
  )
  AUCell::getAUC(auc)
}


#' @rdname BenchGeneScoreSeuratMethod
#' @keywords internal
BenchGeneScoreSeurat.UCell <- function(
  obj,
  pathway,
  assay = NULL,
  layer = "data",
  ..., # pass to UCell
  verbose = getFuncOption("verbose")
) {
  rlang::check_installed("UCell", reason = "to compute UCell scores")
  expr <- SeuratObject::LayerData(object = obj, assay = assay, layer = layer)

  dots <- list(...)

  scores <- do.call(
    UCell::ScoreSignatures_UCell,
    modify_list(
      # `name = ""` keeps gene-set names untouched, matching ssGSEA and AUCell
      list(matrix = expr, features = pathway, name = ""),
      dots
    )
  )

  # UCell returns cells in rows and gene sets in columns
  t(as.matrix(scores))
}


#' @rdname BenchGeneScoreSeuratMethod
#' @keywords internal
BenchGeneScoreSeurat.AddModuleScore <- function(
  obj,
  pathway,
  assay = NULL,
  layer = "data",
  ..., # pass to Seurat::AddModuleScore
  verbose = getFuncOption("verbose")
) {
  rlang::check_installed("Seurat", reason = "to compute AddModuleScore scores")
  assay <- assay %||% SeuratObject::DefaultAssay(obj)
  expr <- SeuratObject::LayerData(object = obj, assay = assay, layer = layer)

  # * AddModuleScore fails on gene sets without any gene in the assay
  pathway <- lapply(X = pathway, FUN = function(genes) genes[genes %chin% rownames(expr)])
  pathway <- pathway[lengths(pathway) > 0L]
  if (length(pathway) == 0L) {
    Abort(
      "None of the gene sets in {.arg pathway} is present in assay {.val {assay}}.",
      "Check the gene names, or the {.arg assay}/{.arg layer} used."
    )
  }

  dots <- list(...)

  args <- modify_list(
    list(
      object = obj,
      features = pathway,
      assay = assay,
      # * internal column prefix, the meta.data of `obj` is left untouched
      name = ".BenchGeneScore",
      slot = layer
    ),
    dots
  )

  scored <- do.call(Seurat::AddModuleScore, args)

  # AddModuleScore stores the score of the i-th gene set in <name><i>
  columns <- paste0(args$name, seq_along(pathway))
  scores <- t(as.matrix(scored[[]][, columns, drop = FALSE]))
  rownames(scores) <- names(pathway)
  scores
}
