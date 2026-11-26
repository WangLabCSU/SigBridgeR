# R Heap High-Water Mark

Returns the peak size of the R heap in megabytes, taken from the
`max used` figures of [`gc()`](https://rdrr.io/r/base/gc.html). A
`gc(reset = TRUE)` before a measurement resets this high-water mark, so
later calls report the peak reached during that measurement rather than
since the process started.

## Usage

``` r
bench_heap_used_mb()
```

## Value

A single `numeric` value: peak R heap usage in MB.
