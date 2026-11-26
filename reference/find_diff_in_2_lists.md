# Find the First Difference Between Two Lists

Recursively compares two lists, matching by name or by position.

## Usage

``` r
find_diff_in_2_lists(x, y)
```

## Arguments

- x, y:

  Two lists to compare.

## Value

`TRUE` (length 1) when they are identical; otherwise a list of three
elements: `name` (path of the first difference), `x` and `y` (the
differing objects, `NULL` when missing on one side).
