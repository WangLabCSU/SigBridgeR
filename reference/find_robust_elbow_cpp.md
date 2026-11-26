# Compute Elbow Heuristics from PCA Embeddings

Compute Elbow Heuristics from PCA Embeddings

## Usage

``` r
find_robust_elbow_cpp(initialized_matrix, ndims)
```

## Arguments

- initialized_matrix:

  A beachmat-compatible matrix whose first `ndims` columns hold PCA
  coordinates.

- ndims:

  `integer`. Number of leading dimensions to consider.

## Value

A list of 11 elements: numeric vectors `stdev`, `variance`,
`pct_variance` and `cumulative_variance` (each of length `ndims`), and
integer scalars `cumulative_90`, `cumulative_80`, `above_mean`,
`above_2sd`, `second_deriv`, `distance` and `final_dims` (the selected
number of dimensions).
