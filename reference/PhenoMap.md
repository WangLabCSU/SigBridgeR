# Map Phenotype Values Using Conditional Rules

Efficiently transforms phenotype values based on user-defined
conditional rules. Similar to
[`dplyr::case_when`](https://dplyr.tidyverse.org/reference/case-and-replace-when.html)
but optimized for performance on large datasets.

## Usage

``` r
PhenoMap(data, ..., .default = NA)
```

## Arguments

- data:

  A vector, data frame, or data.table containing the phenotype data to
  transform.

- ...:

  Named or unnamed formulas of the form `condition ~ value`. Conditions
  are evaluated sequentially; the first matching condition determines
  the output value. Example: `col > 10 ~ "High", col <= 10 ~ "Low"`

- .default:

  Value to use when no conditions match. Default: `NA`.

## Value

A data.table with the transformed column. If input is a vector, returns
a vector.

## See also

[`fcase`](https://rdrr.io/pkg/data.table/man/fcase.html),
[`case_when`](https://dplyr.tidyverse.org/reference/case-and-replace-when.html)

Other input_preprocess:
[`BulkPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/BulkPreProcess.md),
[`PhenoPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/PhenoPreProcess.md),
[`SCPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcess.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Discretize a continuous phenotype
scores <- rnorm(100, mean = 50, sd = 10)
scores <- PhenoMap(scores, scores > 60 ~ "High", scores > 40 ~ "Medium", .default = "Low")

# Example 2: Transform a column in a data frame
df <- data.frame(
  age = c(25, 35, 45, 55, 65),
  gender = c("M", "F", "M", "F", "M")
)
df <- PhenoMap(df, age < 30 ~ "Young", age < 50 ~ "Middle", .default = "Senior")

# Example 3: Multiple conditions with complex logic
df <- data.frame(value = c(5, 15, 25, 35, 45))
df <- PhenoMap(
  df,
  value < 10 ~ "Very Low",
  value < 20 ~ "Low",
  value < 30 ~ "Medium",
  value < 40 ~ "High",
  .default = "Very High"
)
} # }
```
