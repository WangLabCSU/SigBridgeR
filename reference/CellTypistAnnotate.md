# Annotate Cell Types Using CellTypist (Python Backend)

Performs automated cell type annotation via the `celltypist` Python
package using `reticulate` integration. Accepts a `Seurat` object as
input and returns it enriched with CellTypist prediction results as
metadata columns.

Requires a Python environment with `celltypist` installed (see
<https://github.com/Teichlab/celltypist>). The function automatically
attempts to locate a suitable Python interpreter, but users may specify
a custom path via the `python` argument.

## Usage

``` r
CellTypistAnnotate(
  sc,
  model = NULL,
  download = TRUE,
  conda = NULL,
  python = NULL,
  venv_locations = NULL,
  force_update = TRUE,
  verbose = getFuncOption("verbose"),
  celltypist_tools = system.file("python/73-CellTypistAnnotate.py", package =
    "SigBridgeR"),
  ...
)
```

## Arguments

- sc:

  A `Seurat` object containing single-cell RNA-seq data. Must have a
  valid assay with gene expression matrix.

- model:

  Character. CellTypist model specification. One of:

  - Model name (e.g., `"Immune_All_Low"`): Loads the specified built-in
    model.

  - Path to a local `.pkl` model file.

- download:

  Logical. Whether to automatically download the model first. Default:
  `TRUE`.

- conda:

  Character. Conda environment name. Ignore `python` if specified.

- python:

  Character. Path to Python executable with `celltypist` installed. If
  `NULL` (default), auto-detected via
  [`ListPyEnv()`](https://wanglabcsu.github.io/sigbridger/reference/ListPyEnv.md).
  Must point to a valid Python binary.

- venv_locations:

  Character. Path to parent dirtectory storing Python virtual
  environment.

- force_update:

  Logical. Whether to force update the model file. `download` must be
  `TRUE` before this option is effective.

- verbose:

  Logical. Whether to print progress messages during annotation.
  Default: inherits from package option
  `getOption("SigBridgeR.verbose")`.

- celltypist_tools:

  Character. Path to the internal Python bridge script. Default:
  internal package resource
  (`system.file("python/73-CellTypistAnnotate.py", package = "SigBridgeR")`).
  Typically should not be modified by users.

- ...:

  Additional arguments passed to CellTypist's `annotate()` function via
  Python, such as:

  - `majority_voting`: Logical. Whether to refine predicted labels by
    running majority voting classifier after over-clustering. (Default:
    `FALSE`)

  - `mode`: Character. Prediction mode (`"best match"` or
    `"prob match"`). For `"best match"`, selects cell type with largest
    score; `"prob match"` enables multi-label classification. (Default:
    `"best match"`)

  - `p_thres`: Numeric. Probability threshold for multi-label
    classification in `"prob match"` mode. Ignored if
    `mode = "best match"`. (Default: 0.5)

  - `transpose_input`: Logical. Whether to transpose input matrix. Set
    to `TRUE` if filename is in gene-by-cell format. (Default: `FALSE`)

  - `gene_file`: Character. Path to file with genes (one per line)
    corresponding to rows in mtx file. Ignored if input is not in mtx
    format.

  - `cell_file`: Character. Path to file with cells (one per line)
    corresponding to columns in mtx file. Ignored if input is not in mtx
    format.

  - `over_clustering`: Character or vector. Over-clustering
    specification: (1) path to plain file with one cluster ID per
    line; (2) metadata column name in AnnData; (3) vector/array of
    cluster assignments; or (4) omitted for heuristic approach. Ignored
    if `majority_voting = FALSE`.

  - `use_GPU`: Logical. Whether to use GPU acceleration via
    rapids-singlecell for over-clustering. Only relevant when
    `majority_voting = TRUE`. (Default: `FALSE`)

  - `min_prop`: Numeric. Minimum proportion of dominant cell type
    required to name a subcluster. Subclusters below threshold are
    labeled 'Heterogeneous'. Ignored if `majority_voting = FALSE`.
    (Default: 0)

## Value

The input `Seurat` object with some metadata columns added, Column names
may vary slightly depending on CellTypist version and options used.
Usually cell type labels will be added to `meta.data`, and scoring
matrix will be add to `misc$celltypist`

## Requirements

- R packages: `reticulate`, `AnnDataR`

- Python packages: `celltypist`, `scanpy`, `anndata`

- A working Python environment discoverable by `reticulate`

## See also

Other Single_Cell_Annotation_Method:
[`RegisterAnnoMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterAnnoMethod.md),
[`SCAnnotateStrategy`](https://wanglabcsu.github.io/sigbridger/reference/SCAnnotateStrategy.md),
[`SingleRAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/SingleRAnnotate.md),
[`mLLMCelltypeAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/mLLMCelltypeAnnotate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Use a specific immune model with majority voting
annotated <- CellTypistAnnotate(
  seurat_obj,
  model = "Immune_All_Low",
  majority_voting = TRUE
)

# Specify custom Python environment
annotated <- CellTypistAnnotate(
  seurat_obj,
  python = "/path/to/miniconda3/envs/celltypist/bin/python"
)
} # }
```
