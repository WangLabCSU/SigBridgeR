# Normalize and Scale a Matrix Column-wise

Applies the norm and scale transformations to every column of `X`.

## Usage

``` r
normalizeScale_ptr(X)
```

## Arguments

- X:

  Numeric matrix with features in rows and samples in columns.

## Value

Numeric matrix of dimension `ncol(X) x nrow(X)`, i.e. transposed
relative to `X`, matching the original R implementation.
