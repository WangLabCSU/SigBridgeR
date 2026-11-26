# Standardize a Numeric Vector

Computes `(x - mean(x, na.rm = TRUE)) / (sd(x, na.rm = TRUE) + 1e-3)`.

## Usage

``` r
normFunc_ptr(x)
```

## Arguments

- x:

  Numeric vector.

## Value

Numeric vector of the same length as `x`, with missing entries
preserved.
