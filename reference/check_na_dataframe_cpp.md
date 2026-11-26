# Scan a data.frame for Missing Values

Scan a data.frame for Missing Values

## Usage

``` r
check_na_dataframe_cpp(df, nr)
```

## Arguments

- df:

  A data.frame whose columns are scanned for NA/NaN.

- nr:

  `integer`. Expected number of rows; every column must have this
  length.

## Value

A list with `count` (number of missing values) and `row`/`col` (1-based
integer vectors of their positions).
