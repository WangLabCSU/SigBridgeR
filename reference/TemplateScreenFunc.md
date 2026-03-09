# Generate a Template Screening Function with Optional Roxygen2 Documentation

Creates a scaffold for a new screening function that conforms to the
expected interface for bulk-scRNA integration workflows. The template
includes standard parameters (e.g., `matched_bulk`, `sc_data`,
`phenotype`) and a return structure compatible with downstream
validation (e.g., `ValidateScreenFunc`). Optionally includes roxygen2
documentation and supports writing to a new or existing `.R` file.

When `filename = "current"`, the function attempts to append to the
active script in RStudio or Positron (requires `rstudioapi`).

## Usage

``` r
TemplateScreenFunc(
  filename = "my_screen.R",
  func_name = "my_screen",
  append = TRUE,
  documentation = FALSE,
  open = TRUE,
  ...
)
```

## Arguments

- filename:

  Character. Path to the output `.R` file. If `"current"`
  (case-insensitive), uses the currently open script in
  RStudio/Positron. Must end in `.R` if specified explicitly. Default:
  `"my_screen.R"`.

- func_name:

  Character. Name of the generated function. Must be a valid R
  identifier (starts with letter, followed by letters, digits, `_`, or
  `.`). Default: `"my_screen"`.

- append:

  Logical. If `TRUE` (default), appends to an existing file; otherwise,
  overwrites it. If overwriting a non-empty file, prompts for user
  confirmation.

- documentation:

  Logical. Whether to include a full roxygen2 documentation block above
  the function. Default: `FALSE`. If `TRUE`, users should run
  `devtools::document()` afterward.

- open:

  Logical. Whether to open the file in the RStudio/Positron editor after
  writing (if available). Default: `TRUE`.

- ...:

  Additional arguments (currently unused, reserved for future
  extension).

## Value

Invisibly returns the path to the written file (`filename`).

## See also

Other Add_Screen_method:
[`RegisterScreenMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterScreenMethod.md),
[`ScreenStrategy`](https://wanglabcsu.github.io/sigbridger/reference/ScreenStrategy.md),
[`ValidateScreenFunc()`](https://wanglabcsu.github.io/sigbridger/reference/ValidateScreenFunc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a new screen function with documentation
TemplateScreenFunc(
  filename = "custom_screens.R",
  func_name = "MyCoxScreen",
  documentation = TRUE
)

# Append to current script in RStudio
TemplateScreenFunc(filename = "current", func_name = "QuickFilter")
} # }
```
