# Get Names for Dot Variables

Extract or infer names from `...` arguments.

## Usage

``` r
get_names_4_ids(..., .quoses = NULL)
```

## Arguments

- ...:

  Arguments to extract names from.

- .quoses:

  predefined

## Value

Character vector of names.

## Examples

``` r
if (FALSE) { # \dontrun{
get_names_4_ids(a = 1, b = 2)
# [1] "a" "b"
c <- 3
get_names_4_ids(a = 1, b = 2, c)
# [1] "a" "b" "c"
get_names_4_ids()
# character(0)

mat1 <- matrix(1:2)
mat2 <- matrix(3:4)
get_names_4_ids(mat1, mat2)
# [1] "mat1" "mat2"
get_names_4_ids(mat1, c = mat2)
# [1] "mat1" "c"
} # }
```
