#' @title Unified Interface for Single-Cell Annotation Methods
#' @description
#' A dispatch wrapper that routes cell type annotation requests to registered methods
#' stored in \code{SCAnnotateStrategy}. Provides a consistent API for applying
#' different annotation strategies (e.g., \code{CellTypist}, \code{SingleR},
#' \code{mLLMCelltype}) to a \code{Seurat} object without needing to call each
#' method directly.
#'
#' @param sc A \code{Seurat} object containing single-cell RNA-seq data with
#'   preprocessed expression matrix and (optionally) pre-computed clusters.
#' @param method Character. Annotation method to use. Must be one of the keys
#'   registered in \code{SCAnnotateStrategy} (e.g., \code{"CellTypist"},
#'   \code{"SingleR"}, \code{"mLLMCelltype"}). Partial matching is supported.
#'   Use \code{names(SCAnnotateStrategy)} to list available methods.
#' @param ... Method-specific parameters passed to the underlying annotation function.
#'   Examples:
#'   \itemize{
#'     \item For \code{"SingleR"}: \code{ref}, \code{labels}, \code{de.n}, etc.
#'     \item For \code{"CellTypist"}: \code{model}, \code{download}, \code{python}, etc.
#'     \item For \code{"mLLMCelltype"}: \code{tissue_name}, \code{models}, \code{api_keys}, etc.
#'   }
#'   See individual method documentation for details.
#'
#' @return The input \code{Seurat} object enriched with method-specific annotation
#'   metadata columns. Column names and content depend on the chosen method:
#'   \describe{
#'     \item{\code{SingleR}}{\code{SingleR_labels}, \code{SingleR_delta_next}, \code{SingleR_pruned_labels}}
#'     \item{\code{CellTypist}}{\code{celltypist_prediction}, \code{celltypist_probability}, etc.}
#'     \item{\code{mLLMCelltype}}{\code{mllmcelltype_cell_type}, \code{mllmcelltype_consensus_proportion}, \code{mllmcelltype_entropy}}
#'   }
#'
#' @section Extensibility:
#'   New annotation methods can be registered via \code{\link{RegisterAnnoMethod}}
#'   or the unified \code{\link{Register}} interface. Once registered, they become
#'   immediately available through this dispatcher.
#'
#' @examples
#' \dontrun{
#' # List available annotation methods
#' names(SCAnnotateStrategy)
#'
#' # SingleR annotation with HPCA reference
#' annotated <- SCAnnotate(
#'   sc = seurat_obj,
#'   method = "SingleR",
#'   ref = "HPCA"
#' )
#'
#' # CellTypist annotation with custom model
#' annotated <- SCAnnotate(
#'   sc = seurat_obj,
#'   method = "CellTypist",
#'   model = "Immune_All_Low",
#'   download = TRUE
#' )
#'
#' # mLLMCelltype consensus annotation
#' annotated <- SCAnnotate(
#'   sc = seurat_obj,
#'   method = "mLLMCelltype",
#'   tissue_name = "Lung Tumor",
#'   models = c("gpt-4o", "claude-3-5-sonnet-20240620"),
#'   api_keys = list(
#'     openai = Sys.getenv("OPENAI_API_KEY"),
#'     anthropic = Sys.getenv("ANTHROPIC_API_KEY")
#'   )
#' )
#'
#' # Partial matching works
#' annotated <- SCAnnotate(sc = seurat_obj, method = "single")  # -> "SingleR"
#' }
#' @export
#' @family single_cell_preprocess
SCAnnotate <- function(sc, ...) {
  UseMethod("SCAnnotate")
}

SCAnnotate.default <- function(sc, ...) {
  cls_sc <- class(sc)
  Abort(
    "Unsupported class of sc",
    "Expected a {.cls Seurat}, but got a {.cls {cls_sc}}"
  )
}

SCAnnotate.Seurat <- function(
  sc,
  method = c("CellTypist", "SingleR", "mLLMCelltype"),
  ...
) {
  dots <- list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose")
  seed <- dots$verbose %||% getFuncOption("seed")

  set.seed(seed)

  if (length(method) > 1) {
    method <- SigBridgeRUtils::MatchFunc2Args(
      dots,
      mLLMCelltypeAnnotate,
      SingleRAnnotate,
      CellTypistAnnotate,
      name_only = TRUE
    )
    if (length(method) != 1) {
      Abort(
        "x" = "Cannot auto-find a suitable method, please specify a method"
      )
    } else if (verbose) {
      cli::cli_alert_info(
        "Multiple methods detected, auto-using {.pkg {method}}"
      )
    }

    anno_func <- switch(
      method,
      "mLLMCelltypeAnnotate" = mLLMCelltypeAnnotate,
      "SingleRAnnotate" = SingleRAnnotate,
      "CellTypistAnnotate" = CellTypistAnnotate
    )
  } else {
    method <- SigBridgeRUtils::MatchArg(
      method,
      names(SCAnnotateStrategy),
      NULL
    ) # must chosen

    anno_func <- SCAnnotateStrategy[[method]]$executor
  }

  if (method %chin% c("mLLMCelltype", "mLLMCelltypeAnnotate")) {
    anno_func(
      sc,
      ...
    )
  } else if (method %chin% c("SingleR", "SingleRAnnotate")) {
    anno_func(
      sc,
      ref = dots$ref %||% "HPCA",
      ...
    )
  } else if (method %chin% c("CellTypist", "CellTypistAnnotate")) {
    dots <- list2(...)

    if (is.null(dots$conda) && is.null(dots$python)) {
      existing_envs <- ListPyEnv(verbose = FALSE)
      if (!"r-reticulate-celltypist" %chin% existing_envs$name) {
        choice <- utils::askYesNo(
          "Create a new conda environment for CellTypist?"
        )

        if (!isTRUE(choice)) {
          Abort(
            "Aborted. Please specify a conda environment or python interpreter"
          )
        }

        default_args <- list(
          env_type = "conda",
          env_name = "r-reticulate-celltypist",
          method = c("environment"),
          env_file = system.file(
            "conda/celltypist_environment.yml",
            package = "SigBridgeR"
          ),
          python_version = "3.9.15",
          packages = c(
            "celltypist" = "any"
          ),
          env.verbose = SigBridgeRUtils::getFuncOption("verbose")
        )

        exec(
          SetupPyEnv,
          utils::modifyList(
            default_args,
            SigBridgeRUtils::FilterArgs4Func(dots, SetupPyEnv)
          )
        )
      } else if (verbose) {
        ts_cli$cli_alert_info(
          "Existing environment {.val r-reticulate-celltypist} found"
        )
      }
      dots$conda <- "r-reticulate-celltypist"
    }

    # * used to filter out relevant arguments, pass to python function celltypist.annotate
    celltypist.annotate <- function(
      #   filename,
      #   model,
      transpose_input,
      gene_file,
      cell_file,
      mode,
      p_thres,
      majority_voting,
      over_clustering,
      use_GPU,
      min_prop
    ) {
      1
    }

    # * default args
    dots$majority_voting <- dots$majority_voting %||% TRUE

    exec(
      anno_func,
      sc = sc,
      !!!SigBridgeRUtils::FilterArgs4Func(dots, anno_func),
      !!!SigBridgeRUtils::FilterArgs4Func(dots, celltypist.annotate)
    )
  } else {
    cli::cli_abort(c("x" = "Unsupported method"))
  }
}
