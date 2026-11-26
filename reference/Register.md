# Unified Registration Interface for Strategy Methods

Registers screening, preprocessing, or annotation methods into their
respective strategy registries. `Register()` accepts named arguments
where the name becomes the registry key and the value (an S7 object or
function) is dispatched to the appropriate `RegisterImpl` method.

## Usage

``` r
Register(..., overwrite = FALSE, verbose = getFuncOption("verbose"))
```

## Arguments

- ...:

  Named arguments where each name is the registry key and each value is
  an S7 object
  ([ScreenMethod](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethod.md)
  or `AnnotationMethod`) or a function (for `SCPreProcessStrategy`).

- overwrite:

  `logical`. Whether to overwrite an existing entry with the same name.
  Default: `FALSE`.

- verbose:

  `logical`. Whether to print registration success messages. Default:
  `getFuncOption("verbose")`.

## Value

Invisibly `TRUE` on successful registration.

## Details

`Register()` iterates over its named `...` arguments, calling
`RegisterImpl()` for each. `RegisterImpl` is an S7 generic that
dispatches based on the class of the value:

- [ScreenMethod](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethod.md)
  -\> registered into
  [ScreenStrategy](https://wanglabcsu.github.io/sigbridger/reference/ScreenStrategy.md).

- [AnnotationMethod](https://wanglabcsu.github.io/sigbridger/reference/AnnotationMethod.md)
  -\> registered into
  [SCAnnotateStrategy](https://wanglabcsu.github.io/sigbridger/reference/SCAnnotateStrategy.md).

- `function` -\> validated as a Seurat/SeuratObject function and
  registered into
  [SCPreProcessStrategy](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcessStrategy.md).

If the method name already exists in the target registry, the function
aborts unless `overwrite = TRUE`.

## See also

[ScreenStrategy](https://wanglabcsu.github.io/sigbridger/reference/ScreenStrategy.md),
[SCAnnotateStrategy](https://wanglabcsu.github.io/sigbridger/reference/SCAnnotateStrategy.md),
[SCPreProcessStrategy](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcessStrategy.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Register a ScreenMethod
method <- ScreenMethod(
  method_name = "Scissor",
  executor = DoScissor
)
Register(Scissor = method)

# Register a Seurat preprocessing function
Register(h = Seurat::RunHarmony)

# Register with overwrite
Register(Scissor = method, overwrite = TRUE)
} # }
```
