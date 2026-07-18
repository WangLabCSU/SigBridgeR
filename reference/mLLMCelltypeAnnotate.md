# Annotate Cell Types Using Multi-LLM Consensus Approach

Performs automated cell type annotation for Seurat clusters by
leveraging multiple large language models (LLMs) to interpret
cluster-specific marker genes. The function supports both single-model
annotation and multi-model consensus generation, with uncertainty
quantification via agreement metrics.

Workflow:

1.  Identifies top marker genes per cluster (via
    [`Seurat::FindAllMarkers()`](https://satijalab.org/seurat/reference/FindAllMarkers.html)
    or user-provided markers).

2.  Queries specified LLMs with marker gene lists and tissue context.

3.  For multiple models: computes consensus annotations and uncertainty
    metrics (consensus proportion, entropy) via
    `mLLMCelltype::interactive_consensus_annotation()`.

4.  For single model: uses direct annotation via
    `mLLMCelltype::annotate_cell_types()`.

5.  Adds results as metadata columns to the Seurat object.

## Usage

``` r
mLLMCelltypeAnnotate(
  sc,
  seurat_obj_markers = NULL,
  tissue_name = "Human Cancer",
  models = c("gpt-5", "claude-sonnet-4-5-20250929", "gemini-3-pro",
    "qwen-max-2025-01-25"),
  api_keys = list(anthropic = "your-anthropic-key", openai = "your-openai-key", gemini =
    "your-google-key", qwen = "your-qwen-key"),
  ...
)
```

## Arguments

- sc:

  A `Seurat` object with pre-computed clusters (stored in `Idents(sc)`
  or `sc$seurat_clusters`).

- seurat_obj_markers:

  Optional pre-computed marker gene table or list (output of
  [`Seurat::FindAllMarkers()`](https://satijalab.org/seurat/reference/FindAllMarkers.html)).
  If `NULL` (default), markers are computed automatically using
  parameters passed via `...`.

- tissue_name:

  Character. Biological context for annotation (e.g., tissue type,
  disease state). Helps LLMs interpret marker genes appropriately.
  Default: `"Human Tumor"`.

- models:

  Character vector of LLM model identifiers. Supported formats:

  - OpenAI: `"gpt-4o"`, `"gpt-4o-mini"`, etc.

  - Anthropic: `"claude-3-5-sonnet-20240620"`, etc.

  - Google: `"gemini-1.5-pro"`, etc.

  - Alibaba: `"qwen-max"`, `"qwen-plus"`, etc.

  Default:
  `c("gpt-4o", "claude-3-5-sonnet-20240620", "gemini-1.5-pro", "qwen-max")`.
  For single-model mode, only the first model is used.

- api_keys:

  Named list of API keys with provider names as keys:

  `openai`

  :   OpenAI API key

  `anthropic`

  :   Anthropic API key

  `gemini`

  :   Google Cloud API key (with Gemini enabled)

  `qwen`

  :   Alibaba DashScope API key

  Example: `list(openai = "sk-...", anthropic = "sk-ant-...")`.
  **Note**: Default placeholder keys (`"your-xxx-key"`) will fail—users
  must supply valid keys.

- ...:

  Additional arguments passed to downstream functions. Parameters are
  routed as follows:

  To `mLLMCelltype::annotate_cell_types()` and `mLLMCelltype::interactive_consensus_annotation()`:

  :   

      `top_gene_count`

      :   Number of top genes to use per cluster (default: 10L).

      `debug`

      :   Logical. If `TRUE`, prints debugging information.

      `base_urls`

      :   Custom API base URLs: single character string (applied
          globally) or named list with provider-specific URLs (e.g.,
          `list(openai = "...", anthropic = "...")`). Useful for
          proxies, enterprise gateways, or testing environments.

      `controversy_threshold`

      :   Consensus proportion threshold (default: 0.7). Clusters below
          this value are flagged as controversial.

      `entropy_threshold`

      :   Entropy threshold for controversial cluster detection
          (default: 1.0).

      `max_discussion_rounds`

      :   Maximum discussion rounds for controversial clusters (default:
          3).

      `consensus_check_model`

      :   Model used for consensus validation.

      `log_dir`

      :   Directory for log storage (default:
          [`tempdir()`](https://rdrr.io/r/base/tempfile.html)).

      `cache_dir`

      :   Directory for cache storage (default:
          [`tempdir()`](https://rdrr.io/r/base/tempfile.html)).

      `use_cache`

      :   Logical. Whether to use cached results (default: `TRUE`).

      `clusters_to_analyze`

      :   Character/numeric vector of cluster IDs to analyze.
          Non-existent IDs trigger warnings.

      `force_rerun`

      :   Logical. If `TRUE`, bypasses cache and forces re-analysis
          (affects discussion phase only). Default: `FALSE`.

## Value

The input `Seurat` object with the following metadata columns added:

- `mllmcelltype_cell_type`:

  Consensus cell type annotation per cell.

- `mllmcelltype_consensus_proportion`:

  (Multi-model only) Proportion of models agreeing on the assigned label
  (range: 0–1). Higher values indicate stronger consensus.

- `mllmcelltype_entropy`:

  (Multi-model only) Shannon entropy of model predictions. Lower values
  indicate higher confidence (less disagreement among models).

**Note**: Uncertainty metrics are *only added in multi-model mode*
(`length(models) > 1`).

## Requirements

- R packages: `mLLMCelltype`, `plyr`, `Seurat`

- Valid API keys for selected LLM providers (costs may apply)

- Internet connectivity for LLM API calls

## See also

`annotate_cell_types`, `interactive_consensus_annotation`

Other Single_Cell_Annotation_Method:
[`CellTypistAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/CellTypistAnnotate.md),
[`RegisterAnnoMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterAnnoMethod.md),
[`SCAnnotateStrategy`](https://wanglabcsu.github.io/sigbridger/reference/SCAnnotateStrategy.md),
[`SingleRAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/SingleRAnnotate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Multi-model consensus annotation
annotated <- mLLMCellTypeAnnotate(
  sc = pbmc_small,
  tissue_name = "Peripheral Blood Mononuclear Cells",
  models = c("gpt-4o", "claude-3-5-sonnet-20240620"),
  api_keys = list(
    openai = Sys.getenv("OPENAI_API_KEY"),
    anthropic = Sys.getenv("ANTHROPIC_API_KEY")
  ),
  top_gene_count = 15,
  min.pct = 0.25,
  logfc.threshold = 0.5
)

# Usr-level markers
markers_list <- list(
  "0" = c("CD3D", "CD3E", "CD2", "IL7R", "LTB"),
  "1" = c("CD14", "LYZ", "CST3", "MS4A7", "FCGR3A")
)
# Example marker data frame
markers_df <- data.frame(
  cluster = c(0, 0, 0, 1, 1, 1),
  gene = c("CD3D", "CD3E", "CD2", "CD14", "LYZ", "CST3"),
  avg_log2FC = c(2.5, 2.3, 2.1, 3.1, 2.8, 2.5),
  p_val_adj = c(0.001, 0.001, 0.002, 0.0001, 0.0002, 0.0005)
)

# Single-model annotation (faster, no consensus metrics)
annotated <- mLLMCellTypeAnnotate(
  sc = pbmc_small,
  models = "gpt-4o-mini",
  api_keys = list(openai = Sys.getenv("OPENAI_API_KEY"))
)

# Inspect results
table(annotated$mllmcelltype_cell_type)
head(annotated$mllmcelltype_consensus_proportion)  # Only present in multi-model mode
} # }
```
