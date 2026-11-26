# List Available CellTypist Models

Lists the CellTypist pre-trained models available locally. By default it
reads all model files under `~/.celltypist/data/models`, the default
CellTypist models directory (overridable via the `CELLTYPIST_FOLDER`
environment variable). Models are downloaded with
[`DownloadCellTypistModels`](https://wanglabcsu.github.io/sigbridger/reference/DownloadCellTypistModels.md)
and used by
[`CellTypistAnnotate`](https://wanglabcsu.github.io/sigbridger/reference/CellTypistAnnotate.md).

## Usage

``` r
ListCellTypistModels(path = NULL, verbose = getFuncOption("verbose"))
```

## Arguments

- path:

  Character. Path to the CellTypist models directory. Default: `NULL`,
  which resolves to `$CELLTYPIST_FOLDER/data/models` (or
  `~/.celltypist/data/models` when the environment variable is unset).

- verbose:

  Logical. Whether to print the model list to the console. Default:
  `getFuncOption("verbose")`.

## Value

A character vector of available model names (without the `.pkl` suffix).
Invisibly returned when `verbose = TRUE`.

## See also

Other Single_Cell_Annotation_Method:
[`CellTypistAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/CellTypistAnnotate.md),
[`DownloadCellTypistModels()`](https://wanglabcsu.github.io/sigbridger/reference/DownloadCellTypistModels.md),
[`RegisterAnnoMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterAnnoMethod.md),
[`SCAnnotateStrategy`](https://wanglabcsu.github.io/sigbridger/reference/SCAnnotateStrategy.md),
[`SingleRAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/SingleRAnnotate.md),
[`mLLMCelltypeAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/mLLMCelltypeAnnotate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# List all locally available CellTypist models
ListCellTypistModels()

# List models in a custom directory, without printing
ListCellTypistModels(path = "~/my_models", verbose = FALSE)
} # }
```
