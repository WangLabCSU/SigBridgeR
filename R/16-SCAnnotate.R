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
SCAnnotate <- function(
  sc,
  method = c("CellTypist", "SingleR", "mLLMCelltype"),
  ...
) {
  method <- SigBridgeRUtils::MatchArg(
    method,
    names(SCAnnotateStrategy),
    NULL
  ) # must chosen

  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose")
  seed <- dots$verbose %||% getFuncOption("seed")

  set.seed(seed)

  if (method == "mLLMCelltype") {
    SCAnnotateStrategy[[method]](
      sc,
      ...
    )
  } else if (method == "SingleR") {
    SCAnnotateStrategy[[method]](
      sc,
      ref = dots$ref %||% "HPCA",
      ...
    )
  } else if (method == "CellTypist") {
    dots <- rlang::list2(...)

    if (is.null(dots$conda) && is.null(dots$python)) {
      existing_envs <- ListPyEnv()
      if (!"r-reticulate-celltypist" %in% existing_envs$name) {
        choice <- utils::askYesNo(
          "Create a new conda environment for CellTypist?"
        )

        if (!isTRUE(choice)) {
          cli::cli_abort(c(
            "x" = "Aborted. Please specify a conda environment or python interpreter"
          ))
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

        rlang::exec(
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
    celltypist.annotate <- \(
      filename,
      model,
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

    rlang::exec(
      SCAnnotateStrategy[[method]],
      sc = sc,
      !!!SigBridgeRUtils::FilterArgs4Func(dots, SCAnnotateStrategy[[method]]),
      !!!SigBridgeRUtils::FilterArgs4Func(dots, celltypist.annotate)
    )
  }
}
