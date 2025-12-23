# A Function to Replace `tidyr::complete()`

Detects if
[`tidyr::complete()`](https://tidyr.tidyverse.org/reference/complete.html)
is available and uses it if so. Otherwise, uses
[`expand.grid()`](https://rdrr.io/r/base/expand.grid.html) and
[`dplyr::left_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html)
to achieve the same result.

## Usage

``` r
complete_counts(data, ..., fill = list(n = 0))
```
