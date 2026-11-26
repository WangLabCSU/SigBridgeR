# Min-Max Normalize a Numeric Vector

Min-Max Normalize a Numeric Vector

## Usage

``` r
normalize_metric(x, invert = FALSE)
```

## Arguments

- x:

  Numeric vector.

- invert:

  `logical`. Whether to reverse the mapping. Default: `FALSE`.

## Value

Numeric vector of the same length as `x`: non-finite entries become
`NA`, and a constant input maps its finite entries to `0.5`.
