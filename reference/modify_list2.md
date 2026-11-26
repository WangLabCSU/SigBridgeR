# Merge Matching Sub-lists of `val` Against `x`

For every element of `val`, if `x` has an element with the same name and
both are lists, the two are merged recursively with
[`modify_list()`](https://wanglabcsu.github.io/sigbridger/reference/modify_list.md);
otherwise the element of `val` is used unchanged.

## Usage

``` r
modify_list2(x, val)
```

## Arguments

- x:

  A list. Provides the base values for matching sub-lists.

- val:

  A list. Provides the new values.

## Value

A list with the same length and names as `val`.

## Details

Note that the result is built from `val` only: its length and names are
those of `val`, and elements of `x` whose names do not appear in `val`
are **not** carried over. This differs from
[`utils::modifyList()`](https://rdrr.io/r/utils/modifyList.html).

- Names are matched exactly (case-sensitive). If `x` has duplicated
  names, the first match is used.

- Elements of `val` that are unnamed, have an `NA` name or an empty name
  are passed through unchanged.

- Unmodified elements are shared with the inputs, not deep-copied.

- Only the `names` attribute of `val` is kept on the result.

## Examples

``` r
if (FALSE) { # \dontrun{
x <- list(
  a = list(p = 1, q = 2),
  b = 1,
  c = list(z = 0)
)
val <- list(
  a = list(q = 20, r = 30), # both lists: merged recursively
  b = list(new = TRUE),     # `x$b` is not a list: replaced by `val$b`
  d = 4                     # not in `x`: kept as is
)

res <- SigBridgeR:::modify_list2(x, val)
str(res)
# res$a: p = 1, q = 20, r = 30
# res$b: list(new = TRUE)
# res$d: 4
# `x$c` is not in `val`, so it is not in the result
names(res)

# unnamed elements are passed through
SigBridgeR:::modify_list2(list(1), list(list(a = 1), 2))

# non-list inputs are rejected
try(SigBridgeR:::modify_list2(1:3, list()))
} # }
```
