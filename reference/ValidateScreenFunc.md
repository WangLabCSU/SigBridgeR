# Validate Custom Screening Function Compliance

Verifies if a user-provided function meets the interface requirements
for the screening pipeline.

## Usage

``` r
ValidateScreenFunc(func, ...)
```

## Arguments

- func:

  Function to validate. Must accept core parameters and return expected
  structure.

- ...:

  For future updates

## Value

`TRUE` if all checks pass, otherwise terminates with diagnostic report

## See also

Other Add_Screen_method:
[`RegisterScreenMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterScreenMethod.md),
[`ScreenStrategy`](https://wanglabcsu.github.io/sigbridger/reference/ScreenStrategy.md),
[`TemplateScreenFunc()`](https://wanglabcsu.github.io/sigbridger/reference/TemplateScreenFunc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Example compliant function
my_screen <- function(
  sc_data,
  matched_bulk,
  phenotype,
  label_type,
  phenotype_class,
  verbose = FALSE,
  ...
) {
  if (verbose) {
    message("Running custom screen")
  }
  list(
    scRNA_data = sc_data, # Must be Seurat object
    results = data.frame(score = runif(10))
  )
}

ValidateScreenFunc(my_screen)
} # }
```
