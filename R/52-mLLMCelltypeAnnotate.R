#' @title Annotate Cell Types Using Multi-LLM Consensus Approach
#' @family Single_Cell_Annotation_Method
#'
#' @description
#' Performs automated cell type annotation for Seurat clusters by leveraging
#' multiple large language models (LLMs) to interpret cluster-specific marker
#' genes. The function supports both single-model annotation and multi-model
#' consensus generation, with uncertainty quantification via agreement metrics.
#'
#' Workflow:
#'
#' 1. Identifies top marker genes per cluster (via `Seurat::FindAllMarkers()`
#'    or user-provided markers).
#' 2. Queries specified LLMs with marker gene lists and tissue context.
#' 3. For multiple models: computes consensus annotations and uncertainty
#'    metrics (consensus proportion, entropy) via
#'    `mLLMCelltype::interactive_consensus_annotation()`.
#' 4. For single model: uses direct annotation via
#'    `mLLMCelltype::annotate_cell_types()`.
#' 5. Adds results as metadata columns to the Seurat object.
#'
#' @param sc A `Seurat` object with pre-computed clusters (stored in
#'   `Idents(sc)` or `sc$seurat_clusters`).
#' @param seurat_obj_markers Optional pre-computed marker gene table (output of
#'   `Seurat::FindAllMarkers()`). If `NULL` (default), markers are computed
#'   automatically using parameters passed via `...`.
#' @param tissue_name Character. Biological context for annotation (e.g.,
#'   tissue type or disease state), which helps LLMs interpret marker genes
#'   appropriately. Default: `NULL`.
#' @param models Character vector of LLM model identifiers. Supported formats:
#'   * OpenAI: `"gpt-4o"`, `"gpt-4o-mini"`, etc.
#'   * Anthropic: `"claude-3-5-sonnet-20240620"`, etc.
#'   * Google: `"gemini-1.5-pro"`, etc.
#'   * Alibaba: `"qwen-max"`, `"qwen-plus"`, etc.
#'
#'   For single-model mode (`length(models) == 1`), only the first model is
#'   used via `mLLMCelltype::annotate_cell_types()`. Must be provided.
#'   Default: `character()`.
#' @param api_keys Named vector of API keys with provider names as keys:
#'
#'   * **`openai`**: OpenAI API key
#'   * **`anthropic`**: Anthropic API key
#'   * **`gemini`**: Google Cloud API key (with Gemini enabled)
#'   * **`qwen`**: Alibaba DashScope API key
#'
#'   Example: `list(openai = "sk-...", anthropic = "sk-ant-...")`.
#'   Must be provided and named; placeholder keys (e.g., `"your-xxx-key"`)
#'   will fail. Default: `list()`.
#' @param top_gene_count Integer. Number of top marker genes per cluster passed
#'   to the LLMs. Default: `10`.
#' @param debug Logical. Whether to enable debug output from
#'   `mLLMCelltype::annotate_cell_types()`. Default: `FALSE`.
#' @param base_urls Character vector or list of custom base URLs for the LLM
#'   API endpoints (e.g., for proxies or local deployments). Default: `NULL`
#'   (use the default endpoints).
#' @param return_reasoning Logical. Whether to return the reasoning text from
#'   the single-model annotation. Default: `FALSE`.
#' @param controversy_threshold Numeric. Consensus proportion below which a
#'   cluster is considered controversial and triggers further discussion.
#'   Default: `0.7`.
#' @param entropy_threshold Numeric. Entropy below which consensus is
#'   considered reached. Default: `1`.
#' @param max_discussion_rounds Integer. Maximum number of discussion rounds in
#'   the consensus process. Default: `3`.
#' @param consensus_check_model Character. Model used to check consensus during
#'   the discussion. Default: `NULL` (use the first model).
#' @param log_dir Character. Directory for saving the consensus process logs.
#'   Default: `"logs"`.
#' @param cache_dir Character. Directory for caching API responses.
#'   Default: `NULL`.
#' @param use_cache Logical. Whether to use cached API responses.
#'   Default: `TRUE`.
#' @param clusters_to_analyze Character or integer vector. Subset of cluster
#'   IDs to analyze; `NULL` analyzes all clusters. Default: `NULL`.
#' @param force_rerun Logical. Whether to ignore cached results and re-run the
#'   annotation. Default: `FALSE`.
#' @param ... Additional arguments routed to downstream functions:
#'   * Arguments matching `Seurat::FindAllMarkers()` formals (e.g., `min.pct`,
#'     `logfc.threshold`, `test.use`) are passed to marker detection when
#'     `seurat_obj_markers = NULL`.
#'   * `verbose`: Logical. Whether to print progress messages. Default:
#'     `getFuncOption("verbose")`.
#'   * `seed`: Integer. Random seed. Default: `getFuncOption("seed")`.
#'
#' @return The input `Seurat` object with the following metadata columns added:
#'
#'   * **`mllmcelltype_cell_type`**: Consensus cell type annotation per cell.
#'   * **`mllmcelltype_consensus_proportion`**: (Multi-model only) Proportion
#'     of models agreeing on the assigned label (range: 0-1). Higher values
#'     indicate stronger consensus.
#'   * **`mllmcelltype_entropy`**: (Multi-model only) Shannon entropy of model
#'     predictions. Lower values indicate higher confidence (less disagreement
#'     among models).
#'
#'   **Note**: Uncertainty metrics are *only added in multi-model mode*
#'   (`length(models) > 1`).
#'
#' @section Requirements:
#'
#' * R packages: `mLLMCelltype`, `Seurat`
#' * Valid API keys for the selected LLM providers (costs may apply)
#' * Internet connectivity for LLM API calls
#'
#' @examples
#' \dontrun{
#' # Multi-model consensus annotation
#' annotated <- mLLMCelltypeAnnotate(
#'   sc = pbmc_small,
#'   tissue_name = "Peripheral Blood Mononuclear Cells",
#'   models = c("gpt-4o", "claude-3-5-sonnet-20240620"),
#'   api_keys = list(
#'     openai = Sys.getenv("OPENAI_API_KEY"),
#'     anthropic = Sys.getenv("ANTHROPIC_API_KEY")
#'   ),
#'   top_gene_count = 15,
#'   min.pct = 0.25,
#'   logfc.threshold = 0.5
#' )
#'
#' # User-provided markers
#' markers_list <- list(
#'   "0" = c("CD3D", "CD3E", "CD2", "IL7R", "LTB"),
#'   "1" = c("CD14", "LYZ", "CST3", "MS4A7", "FCGR3A")
#' )
#' # Example marker data frame
#' markers_df <- data.frame(
#'   cluster = c(0, 0, 0, 1, 1, 1),
#'   gene = c("CD3D", "CD3E", "CD2", "CD14", "LYZ", "CST3"),
#'   avg_log2FC = c(2.5, 2.3, 2.1, 3.1, 2.8, 2.5),
#'   p_val_adj = c(0.001, 0.001, 0.002, 0.0001, 0.0002, 0.0005)
#' )
#'
#' # Single-model annotation (faster, no consensus metrics)
#' annotated <- mLLMCelltypeAnnotate(
#'   sc = pbmc_small,
#'   models = "gpt-4o-mini",
#'   api_keys = list(openai = Sys.getenv("OPENAI_API_KEY"))
#' )
#'
#' # Inspect results
#' table(annotated$mllmcelltype_cell_type)
#' head(annotated$mllmcelltype_consensus_proportion)  # Only present in multi-model mode
#' }
#' @seealso \code{\link[mLLMCelltype]{annotate_cell_types}},
#'          \code{\link[mLLMCelltype]{interactive_consensus_annotation}}
#' @export
mLLMCelltypeAnnotate <- function(
  sc,
  seurat_obj_markers = NULL,
  tissue_name = NULL, # context
  models = vector("character", 0),
  api_keys = list(),
  top_gene_count = 10,
  debug = FALSE,
  base_urls = NULL,
  return_reasoning = FALSE,
  controversy_threshold = 0.7,
  entropy_threshold = 1,
  max_discussion_rounds = 3,
  consensus_check_model = NULL,
  log_dir = "logs",
  cache_dir = NULL,
  use_cache = TRUE,
  clusters_to_analyze = NULL,
  force_rerun = FALSE,
  ...
) {
  check_installed(c("mLLMCelltype", "dplyr"))
  chk::chk_is(sc, "Seurat")
  chk::chk_vector(api_keys)
  chk::chk_vector(models)
  chk::chk_named(api_keys)
  check_model_key(models = models, api_keys = api_keys)

  dots <- list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose")
  seed <- dots$seed %||% getFuncOption("seed")
  set.seed(seed)

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green(
      "[mLLMCelltype] Start annotating cell types"
    ))
  }

  # Find marker genes for each cluster
  if (is.null(seurat_obj_markers)) {
    if (verbose) {
      ts_cli$cli_alert_info("Find marker genes for each clusters")
    }
    seurat_obj_markers <- exec(
      Seurat::FindAllMarkers,
      object = sc,
      !!!SigBridgeRUtils::FilterArgs4Func(dots, Seurat::FindAllMarkers)
    )
  } else if (verbose) {
    ts_cli$cli_alert_info("Use provided marker genes")
  }

  # Using multiple LLM models to run LLMCelltype annotation.
  if (verbose) {
    ts_cli$cli_alert_info(
      "Large language models cell type Annotating"
    )
  }

  if (length(models) == 1L) {
    # One model prediction
    single_model_results <- exec(
      mLLMCelltype::annotate_cell_types,
      input = seurat_obj_markers,
      tissue_name = tissue_name,
      model = models,
      api_key = api_keys[[1L]],
      top_gene_count = top_gene_count,
      debug = debug,
      base_urls = base_urls,
      return_reasoning = return_reasoning,
    )

    sc$mllmcelltype_cell_type <- dplyr::recode_values(
      x = as.character(SeuratObject::Idents(sc)),
      from = as.character(0L:(length(single_model_results) - 1L)),
      to = single_model_results
    )
    if (verbose) {
      ts_cli$cli_alert_success("Annotation finished with single model")
    }

    return(sc)
  }

  # Multiple model prediction
  consensus_results <- exec(
    mLLMCelltype::interactive_consensus_annotation,
    input = seurat_obj_markers,
    tissue_name = tissue_name, # Provide organizational context.
    models = models,
    api_keys = api_keys,
    top_gene_count = top_gene_count,
    controversy_threshold = controversy_threshold,
    entropy_threshold = entropy_threshold,
    max_discussion_rounds = max_discussion_rounds,
    consensus_check_model = consensus_check_model,
    log_dir = log_dir,
    cache_dir = cache_dir,
    use_cache = use_cache,
    base_urls = base_urls,
    clusters_to_analyze = clusters_to_analyze,
    force_rerun = force_rerun
  )

  # Add annotations to the Seurat object.
  ## Obtain cell type annotations from `consensus_results$final_annotations`.
  cluster_to_celltype_map <- consensus_results$final_annotations

  ## Create a new column for cell type identifiers.
  cell_types <- as.character(Seurat::Idents(sc))
  for (cluster_id in names(cluster_to_celltype_map)) {
    cell_types[cell_types == cluster_id] <- cluster_to_celltype_map[[
      cluster_id
    ]]
  }

  ## Add cell type annotations to the Seurat object.
  sc$mllmcelltype_cell_type <- cell_types

  # Add uncertainty metrics.
  ## Extract detailed consensus results containing metrics.
  consensus_details <- consensus_results$initial_results$consensus_results

  ## Create a data frame containing metrics for each cluster.
  uncertainty_metrics <- data.frame(
    cluster_id = names(consensus_details),
    consensus_proportion = sapply(
      consensus_details,
      function(res) res$consensus_proportion
    ),
    entropy = sapply(consensus_details, function(res) res$entropy)
  )

  # Add uncertainty metrics for each cell.
  ## Match each cell with its corresponding clustering metrics using `seurat_clusters`.
  current_clusters <- sc$seurat_clusters
  sc$mllmcelltype_consensus_proportion <- uncertainty_metrics$consensus_proportion[match(
    current_clusters,
    uncertainty_metrics$cluster_id
  )]
  sc$mllmcelltype_entropy <- uncertainty_metrics$entropy[match(
    current_clusters,
    uncertainty_metrics$cluster_id
  )]

  if (verbose) {
    ts_cli$cli_alert_success("Annotation finished with multiple models")
  }

  sc
}

check_model_key <- function(models = vector(), api_keys = vector()) {
  if (length(models) != length(api_keys)) {
    Abort(
      "[mLLMCelltypeAnnotate()] Number of models does not match number of API keys",
      "Length of models: {length(models)}\nLength of api_keys: {length(api_keys)}"
    )
  }

  api_keys <- unlist(api_keys)
  if (length(api_keys) == 0L) {
    Abort("[mLLMCelltypeAnnotate()] API keys not provided")
  }
}
