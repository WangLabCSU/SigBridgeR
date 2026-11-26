# Scan a beachmat Matrix for Missing Values

Scan a beachmat Matrix for Missing Values

## Usage

``` r
check_na_beachmat_matrix_cpp(initialized_matrix, max_store = -1)
```

## Arguments

- initialized_matrix:

  A beachmat-compatible numeric matrix; unstored sparse entries are
  treated as 0.

- max_store:

  `numeric`. Maximum number of positions kept in the result; a negative
  value (default) keeps all. `count` always reports the total.

## Value

A list with `count` (total number of NA/NaN values), `row` and `col`
(1-based integer vectors of the stored positions).
