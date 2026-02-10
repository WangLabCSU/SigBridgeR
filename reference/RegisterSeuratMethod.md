# Register a Seurat Processing Strategy

Dynamically registers one or more new preprocessing strategies into a
strategy environment (i.e., `SCPreProcessStrategy`). Each strategy is
stored under a user-defined key and wraps a Seurat-compatible function
to conform to the standard interface: `function(object, params)`.

Supports both direct function objects and character specifications of
package-qualified function names (e.g., `h = "Seurat::RunHarmony"`).
This enables safe, lazy resolution of external dependencies and
simplifies pipeline configuration via strings.

## Usage

``` r
RegisterSeuratMethod(
  ...,
  overwrite = FALSE,
  registry = SCPreProcessStrategy,
  verbose = getFuncOption("verbose")
)
```

## Arguments

- ...:

  Named arguments where each name is a **strategy key** (e.g., `"h"`)
  and each value is either:

  - A **function**, or

  - A **character string** specifying a function (e.g.,
    `"Seurat::RunHarmony"` or `"stats::prcomp"`).

  The function will be automatically wrapped to accept `object` and
  `params`, making it compatible with
  [`SCPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcess.md)
  pipelines.

- overwrite:

  Logical. If `FALSE` (default), throws an error when attempting to
  replace an existing strategy. Set to `TRUE` to allow updates.

- registry:

  Environment. The target registry to store strategies. Default:
  `SCPreProcessStrategy`.

- verbose:

  Logical. Whether to print registration messages. Default: inherits
  from `getOption("SigBridgeR.verbose")`.

## Value

Invisibly returns `TRUE` on successful registration.

## See also

Other Registering:
[`Register()`](https://wanglabcsu.github.io/sigbridger/reference/Register.md),
[`RegisterAnnoMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterAnnoMethod.md),
[`RegisterScreenMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterScreenMethod.md)

Other single_cell_preprocess:
[`FilterTumorCell()`](https://wanglabcsu.github.io/sigbridger/reference/FilterTumorCell.md),
[`FindRobustElbow()`](https://wanglabcsu.github.io/sigbridger/reference/FindRobustElbow.md),
[`Pattern2Colname()`](https://wanglabcsu.github.io/sigbridger/reference/Pattern2Colname.md),
[`QCPatternDetect()`](https://wanglabcsu.github.io/sigbridger/reference/QCPatternDetect.md),
[`SCAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/SCAnnotate.md),
[`SCIntegrate()`](https://wanglabcsu.github.io/sigbridger/reference/SCIntegrate.md),
[`SCPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcess.md),
[`SCPreProcessStrategy`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcessStrategy.md),
[`compatible_with_3.0.2()`](https://wanglabcsu.github.io/sigbridger/reference/compatible_with_3.0.2.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Register using a package-qualified name
RegisterSeuratMethod(
  h = "Seurat::RunHarmony",
  registry = SCPreProcessStrategy
)

# Attempting to re-register without `overwrite = TRUE` will fail
# RegisterSeuratMethod(h = Seurat::RunHarmony)  # Error!
} # }
```
