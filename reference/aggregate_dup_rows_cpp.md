# Aggregate Rows of a Matrix by Group

Collapses rows sharing the same 1-based group id into one row using the
requested statistic.

## Usage

``` r
aggregate_dup_rows_cpp(initialized_matrix, group_id, n_groups, method)
```

## Arguments

- initialized_matrix:

  A beachmat-compatible numeric matrix.

- group_id:

  Integer vector of 1-based group ids, one per row.

- n_groups:

  `integer`. Number of groups, i.e. rows of the result.

- method:

  Aggregation method: `"max"`, `"sum"`, `"mean"`, `"median"` or
  `"first"`.

## Value

A numeric matrix with `n_groups` rows and `ncol(initialized_matrix)`
columns.
