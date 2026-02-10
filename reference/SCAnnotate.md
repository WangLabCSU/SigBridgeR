# Unified Interface for Single-Cell Annotation Methods

A dispatch wrapper that routes cell type annotation requests to
registered methods stored in `SCAnnotateStrategy`. Provides a consistent
API for applying different annotation strategies (e.g., `CellTypist`,
`SingleR`, `mLLMCelltype`) to a `Seurat` object without needing to call
each method directly.

## Usage

``` r
SCAnnotate(sc, method = c("CellTypist", "SingleR", "mLLMCelltype"), ...)
```

## Arguments

- sc:

  A `Seurat` object containing single-cell RNA-seq data with
  preprocessed expression matrix and (optionally) pre-computed clusters.

- method:

  Character. Annotation method to use. Must be one of the keys
  registered in `SCAnnotateStrategy` (e.g., `"CellTypist"`, `"SingleR"`,
  `"mLLMCelltype"`). Partial matching is supported. Use
  `names(SCAnnotateStrategy)` to list available methods.

- ...:

  Method-specific parameters passed to the underlying annotation
  function. Examples:

  - For `"SingleR"`: `ref`, `labels`, `de.n`, etc.

  - For `"CellTypist"`: `model`, `download`, `python`, etc.

  - For `"mLLMCelltype"`: `tissue_name`, `models`, `api_keys`, etc.

  See individual method documentation for details.

## Value

The input `Seurat` object enriched with method-specific annotation
metadata columns. Column names and content depend on the chosen method:

- `SingleR`:

  `SingleR_labels`, `SingleR_delta_next`, `SingleR_pruned_labels`

- `CellTypist`:

  `celltypist_prediction`, `celltypist_probability`, etc.

- `mLLMCelltype`:

  `mllmcelltype_cell_type`, `mllmcelltype_consensus_proportion`,
  `mllmcelltype_entropy`

## Extensibility

New annotation methods can be registered via
[`RegisterAnnoMethod`](https://wanglabcsu.github.io/sigbridger/reference/RegisterAnnoMethod.md)
or the unified
[`Register`](https://wanglabcsu.github.io/sigbridger/reference/Register.md)
interface. Once registered, they become immediately available through
this dispatcher.

## See also

Other single_cell_preprocess:
[`FilterTumorCell()`](https://wanglabcsu.github.io/sigbridger/reference/FilterTumorCell.md),
[`FindRobustElbow()`](https://wanglabcsu.github.io/sigbridger/reference/FindRobustElbow.md),
[`Pattern2Colname()`](https://wanglabcsu.github.io/sigbridger/reference/Pattern2Colname.md),
[`QCPatternDetect()`](https://wanglabcsu.github.io/sigbridger/reference/QCPatternDetect.md),
[`RegisterSeuratMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterSeuratMethod.md),
[`SCIntegrate()`](https://wanglabcsu.github.io/sigbridger/reference/SCIntegrate.md),
[`SCPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcess.md),
[`SCPreProcessStrategy`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcessStrategy.md),
[`compatible_with_3.0.2()`](https://wanglabcsu.github.io/sigbridger/reference/compatible_with_3.0.2.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# List available annotation methods
names(SCAnnotateStrategy)

# SingleR annotation with HPCA reference
annotated <- SCAnnotate(
  sc = seurat_obj,
  method = "SingleR",
  ref = "HPCA"
)

# CellTypist annotation with custom model
annotated <- SCAnnotate(
  sc = seurat_obj,
  method = "CellTypist",
  model = "Immune_All_Low",
  download = TRUE
)

# mLLMCelltype consensus annotation
annotated <- SCAnnotate(
  sc = seurat_obj,
  method = "mLLMCelltype",
  tissue_name = "Lung Tumor",
  models = c("gpt-4o", "claude-3-5-sonnet-20240620"),
  api_keys = list(
    openai = Sys.getenv("OPENAI_API_KEY"),
    anthropic = Sys.getenv("ANTHROPIC_API_KEY")
  )
)

# Partial matching works
annotated <- SCAnnotate(sc = seurat_obj, method = "single")  # -> "SingleR"
} # }
```
