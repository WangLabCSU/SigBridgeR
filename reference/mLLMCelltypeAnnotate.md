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
    [`mLLMCelltype::interactive_consensus_annotation()`](https://cafferychen777.github.io/mLLMCelltype/reference/interactive_consensus_annotation.html).

4.  For single model: uses direct annotation via
    [`mLLMCelltype::annotate_cell_types()`](https://cafferychen777.github.io/mLLMCelltype/reference/annotate_cell_types.html).

5.  Adds results as metadata columns to the Seurat object.

## Usage

``` r
mLLMCelltypeAnnotate(sc, ...)

# Default S3 method
mLLMCelltypeAnnotate(sc, ...)

# S3 method for class '`SigBridgeR::ScreenMethodResult`'
mLLMCelltypeAnnotate(sc, ...)

# S3 method for class 'SingleCellExperiment'
mLLMCelltypeAnnotate(sc, ...)

# S3 method for class 'Seurat'
mLLMCelltypeAnnotate(
  sc,
  seurat_obj_markers = NULL,
  tissue_name = NULL,
  models = vector("character", 0L),
  api_keys = list(),
  top_gene_count = 10L,
  debug = FALSE,
  base_urls = NULL,
  return_reasoning = FALSE,
  controversy_threshold = 0.7,
  entropy_threshold = 1L,
  max_discussion_rounds = 3L,
  consensus_check_model = NULL,
  log_dir = "logs",
  cache_dir = NULL,
  use_cache = TRUE,
  clusters_to_analyze = NULL,
  force_rerun = FALSE,
  ...
)
```

## Arguments

- sc:

  A single-cell dataset with pre-computed clusters (stored in
  `Idents(sc)` or `sc$seurat_clusters`), one of:

  - A `Seurat` object (recommended);

  - A `SingleCellExperiment` object, converted to `Seurat` via
    [`Seurat::as.Seurat()`](https://satijalab.github.io/seurat-object/reference/as.Seurat.html)
    before annotation (requires `counts`/`logcounts` assays);

  - A
    [ScreenMethodResult](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodResult.md)
    object, whose `scRNA_data` slot is annotated.

- ...:

  Additional arguments routed to downstream functions:

  - Arguments matching
    [`Seurat::FindAllMarkers()`](https://satijalab.org/seurat/reference/FindAllMarkers.html)
    formals (e.g., `min.pct`, `logfc.threshold`, `test.use`) are passed
    to marker detection when `seurat_obj_markers = NULL`.

  - `verbose`: Logical. Whether to print progress messages. Default:
    `getFuncOption("verbose")`.

  - `seed`: Integer. Random seed. Default: `getFuncOption("seed")`.

- seurat_obj_markers:

  Optional pre-computed marker gene table (output of
  [`Seurat::FindAllMarkers()`](https://satijalab.org/seurat/reference/FindAllMarkers.html)).
  If `NULL` (default), markers are computed automatically using
  parameters passed via `...`.

- tissue_name:

  Character. Biological context for annotation (e.g., tissue type or
  disease state), which helps LLMs interpret marker genes appropriately.
  Default: `NULL`.

- models:

  Character vector of LLM model identifiers. Supported formats:

  - OpenAI: `"gpt-4o"`, `"gpt-4o-mini"`, etc.

  - Anthropic: `"claude-3-5-sonnet-20240620"`, etc.

  - Google: `"gemini-1.5-pro"`, etc.

  - Alibaba: `"qwen-max"`, `"qwen-plus"`, etc.

  For single-model mode (`length(models) == 1`), only the first model is
  used via
  [`mLLMCelltype::annotate_cell_types()`](https://cafferychen777.github.io/mLLMCelltype/reference/annotate_cell_types.html).
  Must be provided. Default:
  [`character()`](https://rdrr.io/r/base/character.html).

- api_keys:

  Named vector of API keys with provider names as keys:

  - **`openai`**: OpenAI API key

  - **`anthropic`**: Anthropic API key

  - **`gemini`**: Google Cloud API key (with Gemini enabled)

  - **`qwen`**: Alibaba DashScope API key

  Example: `list(openai = "sk-...", anthropic = "sk-ant-...")`. Must be
  provided and named; placeholder keys (e.g., `"your-xxx-key"`) will
  fail. Default: [`list()`](https://rdrr.io/r/base/list.html).

- top_gene_count:

  Integer. Number of top marker genes per cluster passed to the LLMs.
  Default: `10`.

- debug:

  Logical. Whether to enable debug output from
  [`mLLMCelltype::annotate_cell_types()`](https://cafferychen777.github.io/mLLMCelltype/reference/annotate_cell_types.html).
  Default: `FALSE`.

- base_urls:

  Character vector or list of custom base URLs for the LLM API endpoints
  (e.g., for proxies or local deployments). Default: `NULL` (use the
  default endpoints).

- return_reasoning:

  Logical. Whether to return the reasoning text from the single-model
  annotation. Default: `FALSE`.

- controversy_threshold:

  Numeric. Consensus proportion below which a cluster is considered
  controversial and triggers further discussion. Default: `0.7`.

- entropy_threshold:

  Numeric. Entropy below which consensus is considered reached. Default:
  `1`.

- max_discussion_rounds:

  Integer. Maximum number of discussion rounds in the consensus process.
  Default: `3`.

- consensus_check_model:

  Character. Model used to check consensus during the discussion.
  Default: `NULL` (use the first model).

- log_dir:

  Character. Directory for saving the consensus process logs. Default:
  `"logs"`.

- cache_dir:

  Character. Directory for caching API responses. Default: `NULL`.

- use_cache:

  Logical. Whether to use cached API responses. Default: `TRUE`.

- clusters_to_analyze:

  Character or integer vector. Subset of cluster IDs to analyze; `NULL`
  analyzes all clusters. Default: `NULL`.

- force_rerun:

  Logical. Whether to ignore cached results and re-run the annotation.
  Default: `FALSE`.

## Value

The input `Seurat` object with the following metadata columns added:

- **`mllmcelltype_cell_type`**: Consensus cell type annotation per cell.

- **`mllmcelltype_consensus_proportion`**: (Multi-model only) Proportion
  of models agreeing on the assigned label (range: 0-1). Higher values
  indicate stronger consensus.

- **`mllmcelltype_entropy`**: (Multi-model only) Shannon entropy of
  model predictions. Lower values indicate higher confidence (less
  disagreement among models).

**Note**: Uncertainty metrics are *only added in multi-model mode*
(`length(models) > 1`).

`SingleCellExperiment` and
[ScreenMethodResult](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodResult.md)
inputs are annotated as `Seurat` objects and returned as such.

## Requirements

- R packages: `mLLMCelltype`, `Seurat`

- Valid API keys for the selected LLM providers (costs may apply)

- Internet connectivity for LLM API calls

## See also

[`annotate_cell_types`](https://cafferychen777.github.io/mLLMCelltype/reference/annotate_cell_types.html),
[`interactive_consensus_annotation`](https://cafferychen777.github.io/mLLMCelltype/reference/interactive_consensus_annotation.html)

Other Single_Cell_Annotation_Method:
[`CellTypistAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/CellTypistAnnotate.md),
[`DownloadCellTypistModels()`](https://wanglabcsu.github.io/sigbridger/reference/DownloadCellTypistModels.md),
[`ListCellTypistModels()`](https://wanglabcsu.github.io/sigbridger/reference/ListCellTypistModels.md),
[`RegisterAnnoMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterAnnoMethod.md),
[`SCAnnotateStrategy`](https://wanglabcsu.github.io/sigbridger/reference/SCAnnotateStrategy.md),
[`SingleRAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/SingleRAnnotate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Multi-model consensus annotation
annotated <- mLLMCelltypeAnnotate(
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

# User-provided markers
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
annotated <- mLLMCelltypeAnnotate(
  sc = pbmc_small,
  models = "gpt-4o-mini",
  api_keys = list(openai = Sys.getenv("OPENAI_API_KEY"))
)

# Inspect results
table(annotated$mllmcelltype_cell_type)
head(annotated$mllmcelltype_consensus_proportion)  # Only present in multi-model mode
} # }
```
