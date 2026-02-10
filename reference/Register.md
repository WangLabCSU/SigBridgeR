# Unified Registration Interface for Strategy Methods

A convenience wrapper that dispatches registration requests to the
appropriate strategy-specific registrar based on the target `registry`.
This function provides a single entry point for registering screening,
preprocessing, or annotation methods without needing to call
[`RegisterScreenMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterScreenMethod.md),
[`RegisterSeuratMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterSeuratMethod.md),
or
[`RegisterAnnoMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterAnnoMethod.md)
directly.

Internally routes to:

- `registry = "auto"` → `detect_registry()`

- `registry = "ScreenStrategy"` →
  [`RegisterScreenMethod`](https://wanglabcsu.github.io/sigbridger/reference/RegisterScreenMethod.md)

- `registry = "SCPreProcessStrategy"` →
  [`RegisterSeuratMethod`](https://wanglabcsu.github.io/sigbridger/reference/RegisterSeuratMethod.md)

- `registry = "SCAnnotateStrategy"` →
  [`RegisterAnnoMethod`](https://wanglabcsu.github.io/sigbridger/reference/RegisterAnnoMethod.md)

## Usage

``` r
Register(
  ...,
  registry = c("auto", "ScreenStrategy", "SCPreProcessStrategy", "SCAnnotateStrategy"),
  verbose = getFuncOption("verbose")
)
```

## Arguments

- ...:

  Arguments passed to the underlying registrar. The exact requirements
  depend on the target `registry`:

  `ScreenStrategy`

  :   Named functions with optional `supported_phenotypes`,
      `parameter_mapper`, etc. (see
      [`RegisterScreenMethod`](https://wanglabcsu.github.io/sigbridger/reference/RegisterScreenMethod.md)).

  `SCPreProcessStrategy`

  :   Named functions or character specifications (e.g.,
      `"h" = "Seurat::RunHarmony"`; see
      [`RegisterSeuratMethod`](https://wanglabcsu.github.io/sigbridger/reference/RegisterSeuratMethod.md)).

  `SCAnnotateStrategy`

  :   Named annotation functions (see
      [`RegisterAnnoMethod`](https://wanglabcsu.github.io/sigbridger/reference/RegisterAnnoMethod.md)).

- registry:

  Character. Target strategy environment for registration. Must be one
  of: `"auto"`, `"ScreenStrategy"`, `"SCPreProcessStrategy"`, or
  `"SCAnnotateStrategy"`. Partial matching is supported (e.g.,
  `"screen"` → `"ScreenStrategy"`).

- verbose:

  Logical. Whether to print registration success messages. Default:
  inherits from package option `getOption("SigBridgeRUtils.verbose")`.

## Value

Invisibly returns `TRUE` on successful registration (via the underlying
registrar).

## See also

[`RegisterScreenMethod`](https://wanglabcsu.github.io/sigbridger/reference/RegisterScreenMethod.md),
[`RegisterSeuratMethod`](https://wanglabcsu.github.io/sigbridger/reference/RegisterSeuratMethod.md),
[`RegisterAnnoMethod`](https://wanglabcsu.github.io/sigbridger/reference/RegisterAnnoMethod.md),
[`SCPreProcessStrategy`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcessStrategy.md)
[`SCAnnotateStrategy`](https://wanglabcsu.github.io/sigbridger/reference/SCAnnotateStrategy.md)
[`ScreenStrategy`](https://wanglabcsu.github.io/sigbridger/reference/ScreenStrategy.md)

Other Registering:
[`RegisterAnnoMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterAnnoMethod.md),
[`RegisterScreenMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterScreenMethod.md),
[`RegisterSeuratMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterSeuratMethod.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Register a screening method for binary/survival phenotypes
Register(
  registry = "ScreenStrategy",
  Scissor = DoScissor,
  supported_phenotypes = c("binary", "survival")
)

# Register a preprocessing step (e.g., Harmony integration)
Register(
  registry = "SCPreProcessStrategy",
  h = "Seurat::RunHarmony"
)

# Register an annotation method
Register(
  registry = "SCAnnotateStrategy",
  my_annot = MyCustomAnnotator
)

# auto detects the target registry
Register(my_annot2 = AnotherAnnotator)
} # }
```
