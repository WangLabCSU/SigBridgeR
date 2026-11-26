# Check Whether a Matrix Looks like Counts

Check Whether a Matrix Looks like Counts

## Usage

``` r
is_counts_matrix_cpp(
  initialized_matrix,
  verbose = FALSE,
  integer_tol = 1e-08,
  min_integer_fraction = 0.95
)
```

## Arguments

- initialized_matrix:

  A beachmat-compatible numeric matrix.

- verbose:

  `logical`. Whether to report the reason on failure. Default: `FALSE`.

- integer_tol:

  `numeric`. Tolerance for treating a value as an integer. Default:
  `1e-8`.

- min_integer_fraction:

  `numeric`. Minimum fraction of values that must be close to integers.
  Default: `0.95`.

## Value

A length-1 logical: `TRUE` when the matrix is non-empty, finite,
non-negative, sufficiently integral and every column has a positive
library size.
