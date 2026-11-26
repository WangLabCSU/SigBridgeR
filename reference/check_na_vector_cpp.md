# Scan a Vector for Missing Values

Scan a Vector for Missing Values

## Usage

``` r
check_na_vector_cpp(x)
```

## Arguments

- x:

  An atomic vector or list to scan for NA/NaN.

## Value

A list with `count` (number of missing values) and `positions` (1-based
integer vector of their positions).
