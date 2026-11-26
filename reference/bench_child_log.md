# Collect Child Process Output

Appends the buffered `stdout`/`stderr` of a
[`callr::r_bg()`](https://callr.r-lib.org/reference/r_bg.html) child to
a rolling log. Only the last `n` lines are kept, so a chatty child
cannot grow the parent heap during a long measurement.

## Usage

``` r
bench_child_log(log, process, n = 100L)
```

## Arguments

- log:

  `character`. Lines collected by earlier calls; the returned vector
  keeps at most `n` of the combined old and new lines.

- process:

  A [`callr::r_bg()`](https://callr.r-lib.org/reference/r_bg.html)
  handle whose output and error pipes are read. Pipe errors (e.g. after
  the child has been reaped) are ignored, so already collected lines are
  preserved.

- n:

  `integer`. Maximum number of lines to keep. Default: `100L`.

## Value

A `character` vector of at most `n` lines.
