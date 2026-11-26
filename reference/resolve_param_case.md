# Resolve `...` Names Against Canonical Parameter Names

Matches the names of a captured `...` list against a set of canonical
parameter names, ignoring case, and writes the matched values into the
calling frame. Matched elements are dropped from the returned list, so
the caller can forward what is left as genuine `...`.

## Usage

``` r
resolve_param_case(dots, canonical, env, lower_inner = FALSE)
```

## Arguments

- dots:

  A list, typically `list(...)` from the calling function. Its names are
  compared against `canonical`.

- canonical:

  A character vector of canonical parameter names, typically
  [`rlang::fn_fmls_names()`](https://rlang.r-lib.org/reference/fn_fmls.html).
  The exact spelling of a matched entry is used as the name of the
  variable assigned in `env`.

- env:

  An environment, typically
  [`rlang::current_env()`](https://rlang.r-lib.org/reference/stack.html),
  in which the canonical parameters live. Matched parameters are written
  here; they are not part of the returned list.

- lower_inner:

  A logical flag. When `TRUE`, the inner names of every matched list are
  lowercased as well, so `list(Alpha = list(MaxIter = 1))` becomes
  `list(Alpha = list(maxiter = 1))`.

## Value

A named list holding the elements of `dots` that were not matched, in
their original order. Its length may be zero.

## Details

This is a helper for functions that accept a fixed set of named
sub-argument lists (for example `bench_config`) but want to tolerate
`BenchRounds` as well as `bench_rounds`.

- Matching is case-insensitive on the **ASCII** range only, and is
  otherwise exact: `foo_bar` does not match `foo-bar`. If `canonical`
  contains several names that differ only in case, the first one wins.

- Each matched value must be a plain list, otherwise an error is raised
  naming the offending argument.

- The current value of the canonical parameter is read from `env` (which
  forces its promise) and merged with the supplied value as
  `modify_list(supplied, current)`: fields present in both take the
  value from `env`, fields only in the supplied list are added, nested
  lists are merged recursively, and `NULL` values are kept. The result
  is then assigned back to the canonical name in `env`.

- When several supplied names match the same canonical parameter, they
  are processed in order, so the last one wins.

- Elements of `dots` that are unnamed, or whose name is `NA` or `""`,
  are never matched and are returned unchanged.

- If `dots` has no names at all, it is returned as is without any work.

- Unmatched elements and untouched values are shared with the input, not
  copied.

## See also

[`modify_list()`](https://wanglabcsu.github.io/sigbridger/reference/modify_list.md)

## Examples

``` r
if (FALSE) { # \dontrun{
f <- function(alpha = list(), ...) {
  dots <- resolve_param_case(
    dots = list(...),
    canonical = rlang::fn_fmls_names(),
    env = rlang::current_env(),
    lower_inner = TRUE
  )
  list(alpha = alpha, dots = dots)
}

# `Alpha` matches `alpha`; inner names are lowercased, `extra` is left over
f(Alpha = list(MaxIter = 10), extra = 1) |> str()

# non-list values are rejected
try(f(Alpha = 10))
} # }
```
