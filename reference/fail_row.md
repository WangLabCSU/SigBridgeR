# Build a Failed-Measurement Row

Creates the one-row table returned by
[`bench_mem_measure()`](https://wanglabcsu.github.io/sigbridger/reference/bench_mem_measure.md)
when a child times out or fails, keeping the same columns as a
successful row but with `NA` memory metrics.

## Usage

``` r
fail_row(method_name, round, position, status, error, started)
```

## Arguments

- method_name:

  `character`. Name of the failed screening method.

- round, position:

  `integer`. Round number and position within the round.

- status:

  `character`. Failure status, e.g. `"error"` or `"timeout"`.

- error:

  `character`. Message describing the failure.

- started:

  `POSIXct` time at which the measurement started, used to compute
  `elapsed`.

## Value

A one-row `data.table`.
