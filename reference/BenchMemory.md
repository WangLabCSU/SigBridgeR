# Benchmark Screening Methods: Wall Time and Memory

Runs every screening method listed in a
[BenchMemoryConfig](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfig.md)
on the same data and measures its wall time and memory footprint with
[`bench::mark()`](https://bench.r-lib.org/reference/mark.html).

## Usage

``` r
BenchMemory(obj, ...)

# S3 method for class '`SigBridgeR::BenchMemoryConfig`'
BenchMemory(obj, ..., verbose = getFuncOption("verbose"))

# Default S3 method
BenchMemory(obj, ...)
```

## Arguments

- obj:

  An object of class
  [BenchMemoryConfig](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfig.md)

- ...:

  Additional arguments passed to methods. Currently unused.

- verbose:

  `logical`. Whether to report progress. Default:
  `getFuncOption("verbose")`.

## Value

A `data.table` with one row per screening method per round. Each row
contains the
[`bench::mark()`](https://bench.r-lib.org/reference/mark.html) columns
(`min`, `median`, `itr/sec`, `mem_alloc`, `n_itr`, `n_gc`, `total_time`,
...), with `expression` replaced by the method name, plus `method_name`,
`round`, `position`, `status`, `error`, `mem_peak_mb`, `rss_peak_mb` and
`elapsed`. Successful measurements have `status = "ok"`; failed or
timed-out measurements are kept with `status` set to `"error"` or
`"timeout"` and `NA` memory metrics, while their benchmark columns are
filled with `NA`.

## Details

The benchmark is organised in *rounds*. Within a round every method is
screened once, and the round order is shuffled by default so that
cold-cache, lazy-loading and package-level caching effects are spread
across methods instead of always penalising the first one. One row is
returned per method per round; summaries across rounds are left to the
caller.

Each measurement runs in its own background R process
([`callr::r_bg()`](https://callr.r-lib.org/reference/r_bg.html)) with
its own copies of the inputs, staged in a temporary directory with
`qs2`. A child process therefore cannot leak package-level caches,
loaded Python runtimes or large result objects into the next
measurement. The child is always terminated on exit so it cannot outlive
the call, and `timeout_sec` bounds the wall time of a single method in a
single round.

Three memory numbers are reported because they mean different things:

- `mem_alloc` - bytes allocated during the run (allocation churn, from
  [`bench::mark()`](https://bench.r-lib.org/reference/mark.html)).

- `mem_peak_mb` - peak R heap high-water mark, from `gc(reset = TRUE)`.

- `rss_peak_mb` - peak process resident set size, recorded only when
  `measure_rss = TRUE` (Linux only).

Before the benchmark starts,
[`bench_probe_resources()`](https://wanglabcsu.github.io/sigbridger/reference/bench_probe_resources.md)
collects a machine resource summary. Set `probe_resources = FALSE` in
the `bench_config` controls (see
[`BenchMemoryConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfigSet.md))
to skip this probe, for example when the workspace is large and the
probe's full [`gc()`](https://rdrr.io/r/base/gc.html) would add
noticeable startup time.

## Benchmark controls (`bench_config`)

- `rounds` (`3L`): number of rounds; every method is screened once per
  round.

- `probe_resources` (`TRUE`): whether to call
  [`bench_probe_resources()`](https://wanglabcsu.github.io/sigbridger/reference/bench_probe_resources.md)
  before staging the inputs. Set to `FALSE` to skip the probe and its
  full [`gc()`](https://rdrr.io/r/base/gc.html).

- `min_iterations`, `max_iterations` (`1L`, `1L`): iteration bounds
  passed to
  [`bench::mark()`](https://bench.r-lib.org/reference/mark.html).

- `order` (`"random"`): `"random"` shuffles methods within each round
  using `seed` plus the round number; `"fixed"` keeps the given order.

- `seed` (`getFuncOption("seed")`): base seed for shuffling and for each
  measurement child.

- `memory`, `filter_gc`, `check`, `min_time`, `relative`, `time_unit`:
  forwarded to
  [`bench::mark()`](https://bench.r-lib.org/reference/mark.html).

- `measure_rss` (`FALSE`): also record peak process RSS (Linux only).

- `on_failure` (`"skip"`): `"skip"` records the failure and continues,
  `"stop"` aborts the whole benchmark.

- `set_threads` (`1L`): number of threads the child process is pinned
  to.

- `timeout_sec` (`10800L`): maximum wall-clock seconds for one method in
  one round; use `Inf` to disable the limit.

## See also

[`BenchMemoryConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfig.md),
[`BenchMemoryConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfigSet.md),
[`bench::mark()`](https://bench.r-lib.org/reference/mark.html)

Other Benchmark:
[`BenchAcc()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAcc.md),
[`BenchGSEA()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEA.md),
[`BenchGeneScore()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScore.md)

## Examples

``` r
if (FALSE) { # \dontrun{
config <- BenchMemoryConfig(
  matched_bulk = bulk,
  sc_data = sc,
  phenotype = pheno,
  screen_method = c("scPP", "scPAS", "LP_SGL"),
  method_param = MethodParamSet(list(
    scPP_param = list(Log2FC_cutoff = 0.585),
    scPAS_param = list(permutation_times = 500)
  )),
  bench_config = list(rounds = 3L, timeout_sec = 1800)
)
res <- BenchMemory(config)
res
res[status == "ok", .(median_sec = median(median)), by = method_name]
} # }
```
