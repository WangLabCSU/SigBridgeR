#' @title Annotate Cell Types Using Multi-LLM Consensus Approach
#' @description
#' Performs automated cell type annotation for Seurat clusters by leveraging
#' multiple large language models (LLMs) to interpret cluster-specific marker genes.
#' The function supports both single-model annotation and multi-model consensus
#' generation, with uncertainty quantification via agreement metrics.
#'
#' Workflow:
#' \enumerate{
#'   \item Identifies top marker genes per cluster (via \code{Seurat::FindAllMarkers()}
#'         or user-provided markers).
#'   \item Queries specified LLMs with marker gene lists and tissue context.
#'   \item For multiple models: computes consensus annotations and uncertainty metrics
#'         (consensus proportion, entropy) via \code{mLLMCelltype::interactive_consensus_annotation()}.
#'   \item For single model: uses direct annotation via \code{mLLMCelltype::annotate_cell_types()}.
#'   \item Adds results as metadata columns to the Seurat object.
#' }
#'
#' @param seurat_obj A \code{Seurat} object with pre-computed clusters
#'   (stored in \code{Idents(seurat_obj)} or \code{seurat_obj$seurat_clusters}).
#' @param seurat_obj_markers Optional pre-computed marker gene table or list (output of
#'   \code{Seurat::FindAllMarkers()}). If \code{NULL} (default), markers are computed
#'   automatically using parameters passed via \code{...}.
#' @param tissue_name Character. Biological context for annotation (e.g., tissue type,
#'   disease state). Helps LLMs interpret marker genes appropriately. Default: \code{"Human Tumor"}.
#' @param models Character vector of LLM model identifiers. Supported formats:
#'   \itemize{
#'     \item OpenAI: \code{"gpt-4o"}, \code{"gpt-4o-mini"}, etc.
#'     \item Anthropic: \code{"claude-3-5-sonnet-20240620"}, etc.
#'     \item Google: \code{"gemini-1.5-pro"}, etc.
#'     \item Alibaba: \code{"qwen-max"}, \code{"qwen-plus"}, etc.
#'   }
#'   Default: \code{c("gpt-4o", "claude-3-5-sonnet-20240620", "gemini-1.5-pro", "qwen-max")}.
#'   For single-model mode, only the first model is used.
#' @param api_keys Named list of API keys with provider names as keys:
#'   \describe{
#'     \item{\code{openai}}{OpenAI API key}
#'     \item{\code{anthropic}}{Anthropic API key}
#'     \item{\code{gemini}}{Google Cloud API key (with Gemini enabled)}
#'     \item{\code{qwen}}{Alibaba DashScope API key}
#'   }
#'   Example: \code{list(openai = "sk-...", anthropic = "sk-ant-...")}.
#'   \strong{Note}: Default placeholder keys (\code{"your-xxx-key"}) will fail—users must supply valid keys.
#'
#' @return The input \code{Seurat} object with the following metadata columns added:
#'   \describe{
#'     \item{\code{mllmcelltype_cell_type}}{Consensus cell type annotation per cell.}
#'     \item{\code{mllmcelltype_consensus_proportion}}{(Multi-model only) Proportion of models
#'           agreeing on the assigned label (range: 0–1). Higher values indicate stronger consensus.}
#'     \item{\code{mllmcelltype_entropy}}{(Multi-model only) Shannon entropy of model predictions.
#'           Lower values indicate higher confidence (less disagreement among models).}
#'   }
#'   \strong{Note}: Uncertainty metrics are \emph{only added in multi-model mode}
#'   (\code{length(models) > 1}).
#'
#' @section Requirements:
#'   \itemize{
#'     \item R packages: \code{mLLMCelltype}, \code{plyr}, \code{Seurat}
#'     \item Valid API keys for selected LLM providers (costs may apply)
#'     \item Internet connectivity for LLM API calls
#'   }
#'
#'
#' @examples
#' \dontrun{
#' # Multi-model consensus annotation
#' annotated <- mLLMCellTypeAnnotate(
#'   seurat_obj = pbmc_small,
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
#' # Usr-level markers
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
#' annotated <- mLLMCellTypeAnnotate(
#'   seurat_obj = pbmc_small,
#'   models = "gpt-4o-mini",
#'   api_keys = list(openai = Sys.getenv("OPENAI_API_KEY"))
#' )
#'
#' # Inspect results
#' table(annotated$mllmcelltype_cell_type)
#' head(annotated$mllmcelltype_consensus_proportion)  # Only present in multi-model mode
#' }
#' @keywords Single_Cell_Annotation_Method
#' @seealso \code{\link[mLLMCelltype]{annotate_cell_types}},
#'          \code{\link[mLLMCelltype]{interactive_consensus_annotation}}
#' @export
mLLMCellTypeAnnotate <- function(
  seurat_obj,
  seurat_obj_markers = NULL,
  tissue_name = "Human Tumor",
  models = c(
    "gpt-5",
    "claude-sonnet-4-5-20250929",
    "gemini-3-pro",
    "qwen-max-2025-01-25"
  ),
  api_keys = list(
    anthropic = "your-anthropic-key",
    openai = "your-openai-key",
    gemini = "your-google-key",
    qwen = "your-qwen-key"
  ),
  ...
) {
  rlang::check_installed(c("mLLMCelltype", "plyr"))
  chk::chk_is(seurat_obj, "Seurat")
  chk::chk_list(api_keys)
  chk::chk_vector(models)
  chk::chk_named(api_keys)
  purrr::walk(c(top_gene_count, controversy_threshold), ~ chk::chk_number)
  check_model_key(models = models, api_keys = api_keys)

  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% getFuncOption("verbose")
  seed <- dots$seed %||% getFuncOption("seed")

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green("Start Annotating Cell Types"))
  }

  # Find marker genes for each cluster
  if (is.null(seurat_obj_markers)) {
    if (verbose) {
      ts_cli$cli_alert_info("Find marker genes for each clusters")
    }
    seurat_obj_markers <- rlang::exec(
      Seurat::FindAllMarkers,
      object = seurat_obj,
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

  if (length(models) == 1) {
    # One model prediction
    single_model_results <- rlang::exec(
      mLLMCelltype::annotate_cell_types,
      input = seurat_obj_markers,
      tissue_name = tissue_name,
      model = trimws(tolower(models)),
      api_key = api_keys[[1]],
      !!!SigBridgeRUtils::FilterArgs4Func(
        dots,
        mLLMCelltype::annotate_cell_types
      )
    )

    if (verbose) {
      ts_cli$cli_alert_success("Annotation Finished")
    }

    seurat_obj$mllmcelltype_cell_type <- plyr::mapvalues(
      x = as.character(SeuratObject::Idents(seurat_obj)),
      from = as.character(0:(length(single_model_results) - 1)),
      to = single_model_results
    )

    return(seurat_obj)
  }
  # Multiple model prediction
  consensus_results <- rlang::exec(
    mLLMCelltype::interactive_consensus_annotation,
    input = seurat_obj_markers,
    tissue_name = tissue_name, # Provide organizational context.
    models = trimws(tolower(models)),
    api_keys = api_keys,
    !!!SigBridgeRUtils::FilterArgs4Func(
      dots,
      mLLMCelltype::interactive_consensus_annotation
    )
  )

  # Add annotations to the Seurat object.
  ## Obtain cell type annotations from `consensus_results$final_annotations`.
  cluster_to_celltype_map <- consensus_results$final_annotations

  ## Create a new column for cell type identifiers.
  cell_types <- as.character(Seurat::Idents(seurat_obj))
  for (cluster_id in names(cluster_to_celltype_map)) {
    cell_types[cell_types == cluster_id] <- cluster_to_celltype_map[[
      cluster_id
    ]]
  }

  ## Add cell type annotations to the Seurat object.
  seurat_obj$mllmcelltype_cell_type <- cell_types

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
  current_clusters <- seurat_obj$seurat_clusters
  seurat_obj$mllmcelltype_consensus_proportion <- uncertainty_metrics$consensus_proportion[match(
    current_clusters,
    uncertainty_metrics$cluster_id
  )]
  seurat_obj$mllmcelltype_entropy <- uncertainty_metrics$entropy[match(
    current_clusters,
    uncertainty_metrics$cluster_id
  )]

  if (verbose) {
    ts_cli$cli_alert_info("Annotation Finished")
  }

  seurat_obj
}

#' @keywords internal
check_model_key <- function(models = vector(), api_keys = list()) {
  if (length(models) != length(api_keys)) {
    cli::cli_abort(c(
      "x" = "[mLLMCelltypeAnnotate()] Number of models does not match number of API keys",
      ">" = "Length of models: {length(models)}",
      ">" = "Length of api_keys: {length(api_keys)}"
    ))
  }

  api_keys <- unlist(api_keys)
  if (
    any(
      api_keys %chin%
        c(
          "your-anthropic-key",
          "your-openai-key",
          "your-google-key",
          "your-qwen-key"
        )
    )
  ) {
    cli::cli_abort(c(
      "x" = "[mLLMCelltypeAnnotate()] API keys not provided"
    ))
  }

  if (length(models) > 0) {
    known_models <- tolower(c(
      # OpenAI
      "gpt-5",
      "gpt-4o-mini",
      "gpt-4.1",
      "gpt-4.1-mini",
      "gpt-4.1-nano",
      "gpt-4-turbo",
      "gpt-3.5-turbo",
      "o1",
      "o1-mini",
      "o1-preview",
      "o1-pro",
      # Anthropic
      "claude-sonnet-4-5-20250929",
      "claude-3-5-sonnet-20241022",
      "claude-3-5-haiku-20241022",
      "claude-3-opus-20240229",
      # DeepSeek
      "deepseek-chat",
      "deepseek-reasoner",
      # Google
      "gemini-3-pro",
      "gemini-3-flash",
      "gemini-2.0-flash",
      "gemini-2.0-flash-lite",
      "gemini-1.5-pro",
      "gemini-1.5-flash",
      # Qwen
      "qwen-max-2025-01-25",
      # Stepfun
      "step-2-mini",
      "step-2-16k",
      "step-1-8k",
      # Zhipu
      "glm-4-plus",
      "glm-3-turbo",
      # MiniMax
      "minimax-text-01",
      # Grok
      "grok-3",
      "grok-3-latest",
      "grok-3-fast",
      "grok-3-fast-latest",
      "grok-3-mini",
      "grok-3-mini-latest",
      "grok-3-mini-fast",
      "grok-3-mini-fast-latest",
      # OpenRouter 格式：provider/model-name
      # OpenAI
      "openai/gpt-5",
      "openai/gpt-4o-mini",
      "openai/gpt-4-turbo",
      "openai/gpt-4",
      "openai/gpt-3.5-turbo",
      # Anthropic
      "anthropic/claude-sonnet-4.5",
      "anthropic/claude-haiku-4",
      "anthropic/claude-opus-4.1",
      # Meta
      "meta-llama/llama-3-70b-instruct",
      "meta-llama/llama-3-8b-instruct",
      "meta-llama/llama-2-70b-chat",
      # Google
      "google/gemini-3-pro",
      "google/gemini-3-flash",
      "google/gemini-1.5-pro-latest",
      "google/gemini-1.5-flash",
      # Mistral
      "mistralai/mistral-large",
      "mistralai/mistral-medium",
      "mistralai/mistral-small",
      # 其他
      "microsoft/mai-ds-r1",
      "perplexity/sonar-small-chat",
      "cohere/command-r",
      "deepseek/deepseek-chat",
      "thudm/glm-z1-32b"
    ))

    cleaned_models <- trimws(tolower(models))
    invalid_mask <- !cleaned_models %in% known_models

    if (any(invalid_mask)) {
      invalid_original <- models[invalid_mask]
      cli::cli_abort(c(
        "x" = "[mLLMCelltypeAnnotate()] Invalid model(s) provided",
        "i" = "The following model(s) are not in the allowed list:",
        ">" = "{invalid_original}",
        "i" = "Allowed models include direct calls (e.g., 'gpt-4o-mini') and OpenRouter format (e.g., 'openai/gpt-4o-mini').",
        "i" = "See {.url https://github.com/cafferychen777/mLLMCelltype} for full model list."
      ))
    }
  }
}
