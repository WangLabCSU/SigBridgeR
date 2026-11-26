# Measure One Screening Method

Runs `method_name` once inside
[`bench::mark()`](https://bench.r-lib.org/reference/mark.html) in an
isolated [`callr::r_bg()`](https://callr.r-lib.org/reference/r_bg.html)
child and returns a one-row table with the timing and memory metrics.
The parent supervises the child: its pipes are drained while it runs,
and it is killed when `bench_config$timeout_sec` is exceeded.

## Usage

``` r
bench_mem_measure(
  config,
  method_name,
  round,
  position,
  tmp_dir,
  verbose = FALSE,
  n_cores = 2L
)
```

## Arguments

- config:

  A
  [BenchMemoryConfig](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfig.md)
  object supplying `label_type`, `phenotype_class` and the
  `bench_config` controls (`seed`, `set_threads`, `measure_rss`,
  `on_failure`, `timeout_sec`, and the arguments forwarded to
  [`bench::mark()`](https://bench.r-lib.org/reference/mark.html)).

- method_name:

  `character`. Name of one screening method; its arguments are taken
  from `config@method_param[[paste0(method_name, "_param")]]`.

- round:

  `integer`. Round number, recorded in the result.

- position:

  `integer`. Position of the method within the round.

- tmp_dir:

  `character`. Directory created by \`bench_tmp_dir()“ holding the
  serialised inputs loaded by the child.

- verbose:

  `logical`. Whether to report progress. Default: `FALSE`.

- n_cores:

  `integer`. Number of CPU cores to configure for the child. Used not in
  the child, but in the parent.

## Value

A one-row `data.table`: the
[`bench::mark()`](https://bench.r-lib.org/reference/mark.html) columns,
with `expression` replaced by `method_name`, plus `method_name`,
`round`, `position`, `status`, `error`, `mem_peak_mb`, `rss_peak_mb` and
`elapsed` (total wall-clock seconds of the child run).

## Details

The child sets the configured seed and threads, loads the inputs from
`tmp_dir`, and runs
[`Screen()`](https://wanglabcsu.github.io/sigbridger/reference/Screen.md)
through [`bench::mark()`](https://bench.r-lib.org/reference/mark.html).
It returns the mark result together with the R heap peak
([`bench_heap_used_mb()`](https://wanglabcsu.github.io/sigbridger/reference/bench_heap_used_mb.md))
and, when `measure_rss` is enabled, the process RSS peak
(\`bench_peak_rss_mb()“).

Failure and timeout are handled according to `on_failure`: `"stop"`
aborts, while `"skip"` returns a row from
[`fail_row()`](https://wanglabcsu.github.io/sigbridger/reference/fail_row.md)
with `status` set to `"error"` or `"timeout"`. The child is always
killed and reaped on exit so it cannot outlive the call.
