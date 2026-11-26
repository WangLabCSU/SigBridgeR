# Aggregate Columns of a Matrix by Group

Collapses columns sharing the same 1-based group id into one column
using the requested statistic.

## Usage

``` r
aggregate_dup_cols_cpp(initialized_matrix, group_id, n_groups, method)
```

## Arguments

- initialized_matrix:

  A beachmat-compatible numeric matrix.

- group_id:

  Integer vector of 1-based group ids, one per column.

- n_groups:

  `integer`. Number of groups, i.e. columns of the result.

- method:

  Aggregation method: `"max"`, `"sum"`, `"mean"`, `"median"` or
  `"first"`.

## Value

A numeric matrix with `nrow(initialized_matrix)` rows and `n_groups`
columns.
