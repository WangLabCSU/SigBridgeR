# Inspect Registered Strategy Environments

Lists and optionally inspects the contents of predefined strategy
environments (e.g., `ScreenStrategy`, `SCPreProcessStrategy`,
`AnnotationStrategy`). Returns a summary of registered methods or
variables, with optional detailed introspection of function signatures,
body structure, and closure properties.

This utility aids in debugging, validation, and exploration of
dynamically registered pipelines—particularly useful when developing or
auditing extensible analysis frameworks.

## Usage

``` r
InterceptStrategy(
  strategy = c("ScreenStrategy", "SCPreProcessStrategy", "AnnotationStrategy"),
  detailed = FALSE,
  ...
)
```

## Arguments

- strategy:

  Character. Name of the strategy environment to inspect. Must be one
  of: `"ScreenStrategy"`, `"SCPreProcessStrategy"`, or
  `"AnnotationStrategy"`. Partial matching is supported.

- detailed:

  Logical. If `TRUE`, performs deep inspection of function objects
  (e.g., argument count, presence of `...`, body length, closure
  environment). Default: `FALSE`.

- ...:

  Additional arguments (currently unused, reserved for future
  extension).

## Value

Invisibly returns a `data.frame`:

- If `detailed = FALSE`: one row per entry in the strategy environment,
  with columns depending on content type (e.g., `method_name`, `func`,
  or `var`).

- If `detailed = TRUE`: augmented with introspection columns such as
  `func_n_args`, `func_has_dots`, `func_body_nlines`, `func_is_closure`,
  etc.

## Examples

``` r
if (FALSE) { # \dontrun{
# List all preprocessing steps
InterceptStrategy("SCPreProcessStrategy")

# Get detailed function metadata for screen methods
details <- InterceptStrategy("ScreenStrategy", detailed = TRUE, verbose = FALSE)
head(details[, grep("executor_", names(details))])
} # }
```
