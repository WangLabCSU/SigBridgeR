# Preprocess Counts and Normalize Column-wise

Applies `1.5^log2(x + 1)` element-wise, then the norm and scale
transformations to every column of `X`.

## Usage

``` r
preprocessCounts_ptr(X)
```

## Arguments

- X:

  Numeric matrix with features in rows and samples in columns.

## Value

Numeric matrix of dimension `ncol(X) x nrow(X)`, i.e. transposed
relative to `X`.
