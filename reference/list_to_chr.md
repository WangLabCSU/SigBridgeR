# Convert a List to a Compact Character Representation

Converts an R `list` into a single character string using the format
`"key: value"`. Named elements use their names as keys, while unnamed
elements use their 1-based positions. Nested lists are converted
recursively and enclosed in parentheses. This method is registered as
`list_to_chr`, so it is called automatically when
[`as.character()`](https://rdrr.io/r/base/character.html) is applied to
a list.

## Usage

``` r
list_to_chr(x)
```

## Arguments

- x:

  A list to be converted.

## Value

A character vector of length one.

## Examples

``` r
x <- list(key1 = "val1", key2 = "val2")
as.character(x)
#> [1] "val1" "val2"
# "key1: val1, key2: val2"
x <- list(key1 = list(key2 = list(key3 = "val3")))
as.character(x)
#> [1] "list(key2 = list(key3 = \"val3\"))"
# "key1: (key2: (key3: val3))"
x <- list(
  a = 1,
  b = list(
    c = 2,
    d = list(e = "hello")
  ),
  f = c(1, 2, 3)
)
as.character(x)
#> [1] "1"                                   
#> [2] "list(c = 2, d = list(e = \"hello\"))"
#> [3] "c(1, 2, 3)"                          
# "a: 1, b: (c: 2, d: (e: hello)), f: 1, 2, 3"
x <- list("foo", bar = NA, baz = NULL)
as.character(x)
#> [1] "foo"  "NA"   "NULL"
# "1: foo, bar: NA, baz: NULL"
```
