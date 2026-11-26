# Download CellTypist Models

Downloads CellTypist pre-trained models via the celltypist Python
backend. A single model, multiple models, or all models can be
downloaded. The models are saved to CellTypist's local models directory
and can subsequently be used by
[`CellTypistAnnotate`](https://wanglabcsu.github.io/sigbridger/reference/CellTypistAnnotate.md)
without re-downloading each time.

To see all available models and their descriptions, see
<https://github.com/Teichlab/celltypist>.

## Usage

``` r
DownloadCellTypistModels(model = NULL, force_update = FALSE)
```

## Arguments

- model:

  Character. Name(s) of the CellTypist model(s) to download, e.g.
  `"Immune_All_Low"`. A `.pkl` suffix is appended automatically if
  missing. When `NULL`, all models are downloaded. Default: `NULL`.

- force_update:

  Logical. Whether to force re-downloading the CellTypist model files
  even when they already exist locally. Default: `FALSE`.

## Value

Invisibly `NULL`. The models are downloaded to the CellTypist models
directory (see `models.models_path` in the Python backend).

## Requirements

- A Python environment with `celltypist` installed, discoverable by
  `reticulate`.

## See also

Other Single_Cell_Annotation_Method:
[`CellTypistAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/CellTypistAnnotate.md),
[`ListCellTypistModels()`](https://wanglabcsu.github.io/sigbridger/reference/ListCellTypistModels.md),
[`RegisterAnnoMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterAnnoMethod.md),
[`SCAnnotateStrategy`](https://wanglabcsu.github.io/sigbridger/reference/SCAnnotateStrategy.md),
[`SingleRAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/SingleRAnnotate.md),
[`mLLMCelltypeAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/mLLMCelltypeAnnotate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Download all CellTypist models
DownloadCellTypistModels()

# Download a single model
DownloadCellTypistModels(model = "Immune_All_Low")

# Force re-downloading a specific model
DownloadCellTypistModels(model = "Immune_All_Low", force_update = TRUE)
} # }
```
