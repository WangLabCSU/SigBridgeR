# Annotate Cell Types Using CellTypist (Python Backend)

Performs automated cell type annotation via the `celltypist` Python
package using `reticulate` integration. Accepts a `Seurat` object as
input and returns it enriched with CellTypist prediction results as
metadata columns.

Requires a Python environment with `celltypist` installed (see
<https://github.com/Teichlab/celltypist>). By default the function
locates a suitable Python interpreter automatically via `reticulate`;
use
[`reticulate::use_python()`](https://rstudio.github.io/reticulate/reference/use_python.html)
to pin a specific environment beforehand.

## Usage

``` r
CellTypistAnnotate(sc, ...)

# Default S3 method
CellTypistAnnotate(sc, ...)

# S3 method for class '`SigBridgeR::ScreenMethodResult`'
CellTypistAnnotate(sc, ...)

# S3 method for class 'SingleCellExperiment'
CellTypistAnnotate(sc, ...)

# S3 method for class 'Seurat'
CellTypistAnnotate(
  sc,
  model = NULL,
  download = TRUE,
  force_update = FALSE,
  majority_voting = FALSE,
  mode = c("best match", "prob match"),
  p_thres = 0.5,
  transpose_input = FALSE,
  gene_file = NULL,
  cell_file = NULL,
  over_clustering = NULL,
  use_GPU = FALSE,
  min_prop = 0L,
  verbose = getFuncOption("verbose"),
  assay = "RNA",
  ...
)
```

## Arguments

- sc:

  A single-cell dataset, one of:

  - A `Seurat` object (recommended) with a valid assay and gene
    expression matrix;

  - A `SingleCellExperiment` object, converted to `Seurat` via
    [`Seurat::as.Seurat()`](https://satijalab.github.io/seurat-object/reference/as.Seurat.html)
    before annotation (requires `counts`/`logcounts` assays);

  - A
    [ScreenMethodResult](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodResult.md)
    object, whose `scRNA_data` slot is annotated.

- ...:

  Additional arguments passed through to the CellTypist Python backend.

- model:

  Character. CellTypist model specification. One of:

  - A built-in model name (e.g., `"Immune_All_Low"`); a `.pkl` suffix is
    appended automatically if missing.

  - A path to a local `.pkl` model file.

- download:

  Logical. Whether to automatically download the model first. Default:
  `TRUE`.

- force_update:

  Logical. Whether to force re-downloading the model files when
  `download = TRUE`. Default: `FALSE`.

- majority_voting:

  Logical. Whether to refine predicted labels by running a majority
  voting classifier after over-clustering. Default: `FALSE`.

- mode:

  Character. Prediction mode: `"best match"` selects the cell type with
  the largest score; `"prob match"` enables multi-label classification.
  Default: `"best match"`.

- p_thres:

  Numeric. Probability threshold for multi-label classification in
  `"prob match"` mode. Ignored when `mode = "best match"`. Default:
  `0.5`.

- transpose_input:

  Logical. Whether to transpose the input matrix. Set to `TRUE` if the
  input is in gene-by-cell format. Default: `FALSE`.

- gene_file:

  Character. Path to a file with genes (one per line) corresponding to
  rows of the input matrix. Default: `NULL`.

- cell_file:

  Character. Path to a file with cells (one per line) corresponding to
  columns of the input matrix. Default: `NULL`.

- over_clustering:

  Character or vector. Over-clustering specification used for majority
  voting: (1) path to a plain file with one cluster ID per line; (2)
  metadata column name in AnnData; (3) vector/array of cluster
  assignments; or (4) omitted for the heuristic approach. Ignored when
  `majority_voting = FALSE`. Default: `NULL`.

- use_GPU:

  Logical. Whether to use GPU acceleration (via rapids-singlecell) for
  over-clustering. Only relevant when `majority_voting = TRUE`. Default:
  `FALSE`.

- min_prop:

  Numeric. Minimum proportion of the dominant cell type required to name
  a subcluster; subclusters below the threshold are labeled
  `'Heterogeneous'`. Ignored when `majority_voting = FALSE`. Default:
  `0L`.

- verbose:

  Logical. Whether to print progress messages during annotation.
  Default: `getFuncOption("verbose")`.

- assay:

  which Seurat assay to use

## Value

The input `Seurat` object with prediction results added: cell type
labels are stored as new columns in `meta.data` (prefixed with
`celltypist_`), and the scoring matrices are stored in
`misc$celltypist`. Column names may vary slightly depending on the
CellTypist version and options used. `SingleCellExperiment` and
[ScreenMethodResult](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodResult.md)
inputs are annotated as `Seurat` objects and returned as such.

## Requirements

- R packages: `reticulate`, `AnnDataR`

- Python packages: `celltypist`, `scanpy`, `anndata`

- A working Python environment discoverable by `reticulate`

## See also

Other Single_Cell_Annotation_Method:
[`DownloadCellTypistModels()`](https://wanglabcsu.github.io/sigbridger/reference/DownloadCellTypistModels.md),
[`ListCellTypistModels()`](https://wanglabcsu.github.io/sigbridger/reference/ListCellTypistModels.md),
[`RegisterAnnoMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterAnnoMethod.md),
[`SCAnnotateStrategy`](https://wanglabcsu.github.io/sigbridger/reference/SCAnnotateStrategy.md),
[`SingleRAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/SingleRAnnotate.md),
[`mLLMCelltypeAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/mLLMCelltypeAnnotate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Use a built-in immune model with majority voting
annotated <- CellTypistAnnotate(
  sc = seurat_obj,
  model = "Immune_All_Low",
  majority_voting = TRUE
)

# Probability-based multi-label annotation
annotated <- CellTypistAnnotate(
  sc = seurat_obj,
  model = "Immune_All_Low",
  mode = "prob match",
  p_thres = 0.5
)

# Use a local model file
annotated <- CellTypistAnnotate(
  sc = seurat_obj,
  model = "/path/to/My_Model.pkl"
)
} # }
```
