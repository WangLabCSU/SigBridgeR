# Trace and compute the value of a variable defined inside a function

Recursively traces variable assignments (both from function parameters
and function body) to reconstruct and evaluate the expression that
produces the variable's value. The trace follows the assignment chain
until it reaches literal values or unresolvable symbols.

## Usage

``` r
get_var_value(var_name, func, .env = rlang::current_env())
```

## Arguments

- var_name:

  Character string, the name of the target variable.

- func:

  An R function inside which the variable is defined.

- .env:

  Evaluation environment used for resolving function calls (e.g. `runif`
  from the **stats** package). Defaults to the caller's environment
  ([`rlang::current_env()`](https://rlang.r-lib.org/reference/stack.html)),
  which is typically `.GlobalEnv` and has access to all attached
  packages.

## Value

The computed value of the variable.

## Examples

``` r
f <- function(save_path = "./analysis") {
  save_path_new <- file.path(save_path, "res")
  return(save_path_new)
}
get_var_value("save_path_new", f)   # returns "./analysis/res"
#> [1] "./analysis/res"

g <- function(a = 1, b = 2){
  c <- a * 2 + b * 3
  d = c^2
  e <<- d - 1
  e
}
get_var_value("e", g)   # returns 63
#> [1] 63

h <- function(a = "A", ...){
  a <- 1
  return(a)
  a <- 2
  a
}
get_var_value("a",h)    # returns 1 (dead code after return is ignored)
#> [1] 1

i <- function(x = 2) {
  for (k in 1:3) x <- x * 2
  x
}
get_var_value("x", i)   # returns 16
#> [1] 16

j <- function(x = 2) {
  while (x < 10) x <- x * 2
  x
}
get_var_value("x", j)   # returns 16
#> [1] 16
```
