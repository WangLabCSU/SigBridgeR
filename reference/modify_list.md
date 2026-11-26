# Faster modifyList()

Modifies a possibly nested list recursively by changing a subset of
elements at each level to match a second list.

## Usage

``` r
modify_list(x, val, keep_null = TRUE)
```

## Arguments

- x:

  A named [`list`](https://rdrr.io/r/base/list.html), possibly empty.

- val:

  A named list with components to replace corresponding components in
  `x` or add new components.

- keep_null:

  If TRUE, NULL elements in val become NULL elements in x. Otherwise,
  the corresponding element, if present, is deleted from x.

## Value

A modified version of `x`, with the modifications determined as follows
(here, list elements are identified by their names). Elements in `val`
which are missing from `x` are added to `x`. For elements that are
common to both but are not both lists themselves, the component in `x`
is replaced (or possibly deleted, depending on the value of `keep.null`)
by the one in `val`. For common elements that are in both lists,
`x[[name]]` is replaced by `modifyList(x[[name]], val[[name]])`.

## Author

Deepayan Sarkar <Deepayan.Sarkar@R-project.org>

## Examples

``` r
foo <- list(a = 1, b = list(c = "a", d = FALSE))
bar <- modifyList(foo, list(e = 2, b = list(d = TRUE)))
str(foo)
#> List of 2
#>  $ a: num 1
#>  $ b:List of 2
#>   ..$ c: chr "a"
#>   ..$ d: logi FALSE
str(bar)
#> List of 3
#>  $ a: num 1
#>  $ b:List of 2
#>   ..$ c: chr "a"
#>   ..$ d: logi TRUE
#>  $ e: num 2
```
