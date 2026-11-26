# Test Whether a Target Value Is Over- or Under-Represented

Test Whether a Target Value Is Over- or Under-Represented

## Usage

``` r
IsSkewedDynamic_cpp(x, target = 0, expected_p = 0.8, n_sd = 4L)
```

## Arguments

- x:

  Vector of values to examine.

- target:

  `numeric`. Value to test. Default: `0`.

- expected_p:

  `numeric`. Expected proportion of `target`. Default: `0.8`.

- n_sd:

  `integer`. Number of binomial standard errors allowed. Default: `4L`.

## Value

A length-1 logical: `TRUE` when the observed proportion deviates by more
than `n_sd` standard errors; `NA` for empty input.
