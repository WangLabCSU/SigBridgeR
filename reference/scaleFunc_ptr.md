# Scale a Numeric Vector to 0-1

Computes `(x - min(x)) / (max(x) - min(x) + 1e-3)`; typically applied
after
[`normFunc_ptr()`](https://wanglabcsu.github.io/sigbridger/reference/normFunc_ptr.md).

## Usage

``` r
scaleFunc_ptr(x)
```

## Arguments

- x:

  Numeric vector.

## Value

Numeric vector of the same length as `x`.
