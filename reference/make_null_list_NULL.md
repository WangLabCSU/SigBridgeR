# Replace Empty Lists with NULL Recursively

Replace Empty Lists with NULL Recursively

## Usage

``` r
make_null_list_NULL(x)
```

## Arguments

- x:

  Any R object; lists are processed recursively.

## Value

`x` with every empty list replaced by `NULL`; non-list objects are
returned unchanged and names of non-empty lists are kept.
