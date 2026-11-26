# Convert a List to an Indented Multi-line Character Representation

Converts an R `list` into a single character string in which every
element occupies its own line, using the format `"key: value"`. Named
elements use their names as keys, while unnamed elements use their
1-based positions. Nested lists are expanded recursively, and each
nesting level is indented by two additional spaces.

## Usage

``` r
list_to_chr_formated(x)
```

## Arguments

- x:

  A list to be converted.

## Value

A character vector of length one, with lines separated by `"\n"` (no
trailing newline).

## Examples

``` r
x <- list(
  a = 1,
  b = list(
    c = 2,
    d = list(e = "hello")
  ),
  f = c(1, 2, 3)
)
cat(list_to_chr_formated(x), "\n")
#> a: 1
#> b:
#>   c: 2
#>   d:
#>     e: hello
#> f: 1, 2, 3 
# a: 1
# b:
#   c: 2
#   d:
#     e: hello
# f: 1, 2, 3
```
