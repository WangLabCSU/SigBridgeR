# Register a Custom Screening Method for Phenotype-Driven Analysis

Registers one or more user-defined screening functions into a shared
registry (i.e., `ScreenStrategy`), enabling dynamic dispatch based on
phenotype type (binary, survival, or continuous). Each method is stored
with metadata including supported phenotypes and an optional parameter
remapping function.

This function supports community-driven extension of analysis
pipelines—ideal for integrating novel single cell phenotypic screening
methods while maintaining compatibility with downstream validation and
execution frameworks

## Usage

``` r
RegisterScreenMethod(
  ...,
  supported_phenotypes = c("binary", "survival", "continuous"),
  parameter_mapper = NULL,
  registry = ScreenStrategy,
  verbose = getFuncOption("verbose") %||% TRUE
)
```

## Arguments

- ...:

  Named arguments where each value is a **function** implementing a
  screening method. The name becomes the method identifier (e.g.,
  `Scissor = DoScissor` registers `DoScissor` under the key
  `"Scissor"`). Unnamed arguments are auto-named using their expression
  (e.g., `Scissor` → name `"Scissor"`).

- supported_phenotypes:

  Character vector specifying which phenotype types the method supports.
  Must be a subset of `c("binary", "survival", "continuous")`. Default:
  all three.

- parameter_mapper:

  A function that transforms the input parameter list before passing it
  to the executor. Useful for changing parameters from interface
  function. Receives a named list `params` and must return a modified
  list. Default: `NULL`.

- registry:

  An environment used as a method registry (e.g., `ScreenStrategy`). New
  methods are stored as `registry[["MethodName"]]` = metadata list.
  Default: `ScreenStrategy`.

- verbose:

  Logical. Whether to print a success message upon registration.
  Default: `TRUE`.

## Value

Invisibly returns `TRUE` on successful registration.

## Registry Entry Structure

Each registered method is stored as a list with four elements:

- `method_name`:

  Character: the method key (e.g., "Scissor")

- `executor`:

  Function: the actual screening implementation

- `phenotypes`:

  Character vector: supported phenotype classes

- `mapper`:

  Function: parameter transformation hook

## Examples

``` r
if (FALSE) { # \dontrun{
# Example: Register a custom screen method for survival data
DoMyScreen <- function(...) {
  stop("Not implemented")
}

RegisterScreenMethod(
  # Key = function
  MyCox = DoMyScreen,
  # If `continuous` phenotype is not supported
  supported_phenotypes = c("survival", "binary"),
  parameter_mapper = function(params) {
    # If I perfer to use `quiet` instead of `verbose`
    params$quiet <- !params$verbose
    params
  }
)
} # }
```
