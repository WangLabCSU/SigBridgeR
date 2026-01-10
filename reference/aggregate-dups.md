# Aggregate Rows or Columns with Duplicate Names

These functions collapse duplicated row names (e.g., gene symbols) or
column names (e.g., sample IDs) in matrix-like objects by aggregating
values using configurable methods. They support:

- Rows:

  `AggregateDupRows`: merges rows sharing the same row name.

- Columns:

  `AggregateDupCols`: merges columns sharing the same column name.

- Both:

  `AggregateDups`: convenience wrapper applying row-then-column
  aggregation.

Designed for expression matrices, count tables, or any numeric data
where feature/sample duplication occurs. Handles `matrix`, `data.frame`,
and S4 `Matrix` classes (e.g. `dgCMatrix`) robustly.

Convenience wrapper that first aggregates duplicated rows, then
duplicated columns. Useful for cleaning matrices where both feature and
sample duplication may occur.

## Usage

``` r
AggregateDupRows(
  x,
  method = c("max", "sum", "mean", "median", "first"),
  verbose = TRUE,
  ...
)

AggregateDupCols(
  x,
  method = c("max", "sum", "mean", "median", "first"),
  verbose = TRUE,
  ...
)

AggregateDups(
  x,
  method = c("max", "sum", "mean", "median", "first"),
  row_method = NULL,
  col_method = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  A numeric matrix-like object (see Details).

- method:

  Character scalar. Aggregation method (see Methods below).

- verbose:

  Whether to print messages

- ...:

  No usage

- row_method:

  Aggregation method for rows. Defaults to `method`.

- col_method:

  Aggregation method for columns. Defaults to `method`.

## Value

An aggregated object of the same effective type as `x`, with unique
row/column names.

## Methods

Supported methods (applied column-wise for rows, row-wise for columns):

- `"max"`:

  Maximum value per group (default).

- `"sum"`:

  Sum of values per group.

- `"mean"`:

  Arithmetic mean (uses `na.rm = TRUE`).

- `"median"`:

  Median value.

- `"first"`:

  First occurrence in original order.

## Input Types and Return Types

|              |                                                         |
|--------------|---------------------------------------------------------|
| Input class  | Output class (unless noted)                             |
| `matrix`     | `matrix`                                                |
| `data.frame` | `data.frame`                                            |
| S4 `Matrix`  | `matrix` (dense) — S4 attributes dropped for generality |

Row/column order in output follows *first occurrence* of each unique
name in `rownames(x)` / `colnames(x)`.

## Examples

``` r
# Full deduplication in one step
mat <- matrix(1:16, nrow = 4,
              dimnames = list(c("TP53", "TP53", "BRCA1", "ACTB"),
                            c("S1", "S1", "S2", "S3")))
AggregateDups(mat, method = "sum")
#> ✔ No duplicated column names found.
#>       S1 S1.1 S2 S3
#> TP53   3    3 19 27
#> BRCA1  3    3 11 15
#> ACTB   4    4 12 16
#>       S1 S2 S3
#> TP53   5  7  9
#> BRCA1  3  7 11
#> ACTB   4  8 12
```
