# Register an Annotation Method into the Strategy Registry

Registers one or more user-defined annotation functions into a shared
registry environment (e.g., `SCAnnotateStrategy`). Each method is stored
with minimal metadata (`method_name` and `executor`) to enable dynamic
dispatch in annotation pipelines.

This function supports extensible cell type annotation frameworks—ideal
for integrating novel classifiers (e.g., custom SingleR wrappers,
marker-based assigners, or deep learning models) while maintaining a
uniform execution interface.

## Usage

``` r
RegisterAnnoMethod(
  ...,
  registry = SCAnnotateStrategy,
  overwrite = FALSE,
  verbose = getFuncOption("verbose")
)
```

## Arguments

- ...:

  Named arguments where each value is a function implementing an
  annotation method. The name becomes the method identifier (e.g.,
  `SingleR = SingleRAnnotate` registers under key `"SingleR"`). Unnamed
  arguments are auto-named using their expression via
  [`get_names_4_ids()`](https://wanglabcsu.github.io/sigbridger/reference/get_names_4_ids.md).

- registry:

  Target registry to store annotation methods. Default:
  `SCAnnotateStrategy`.

- overwrite:

  Logical. If `FALSE` (default), throws an error when attempting to
  replace an existing method. Set to `TRUE` to allow updates.

- verbose:

  Logical. Whether to print a success message upon registration.
  Default: inherits from package option
  `getOption("SigBridgeRUtils.verbose")`.

## Value

Invisibly returns `TRUE` on successful registration.

## See also

Other Registering:
[`Register()`](https://wanglabcsu.github.io/sigbridger/reference/Register.md),
[`RegisterScreenMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterScreenMethod.md),
[`RegisterSeuratMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterSeuratMethod.md)

Other Single_Cell_Annotation_Method:
[`CellTypistAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/CellTypistAnnotate.md),
[`SCAnnotateStrategy`](https://wanglabcsu.github.io/sigbridger/reference/SCAnnotateStrategy.md),
[`SingleRAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/SingleRAnnotate.md),
[`mLLMCelltypeAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/mLLMCelltypeAnnotate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Register a custom annotation wrapper
MyAnnotator <- function(sc, ...) {
  # Custom logic returning annotated Seurat/SCE object
  sc
}

RegisterAnnoMethod(
  custom_annot = MyAnnotator,
  registry = SCAnnotateStrategy
)

# Attempting to re-register without `overwrite = TRUE` will fail
# RegisterAnnoMethod(custom_annot = AnotherAnnotator)  # Error!
} # }
```
