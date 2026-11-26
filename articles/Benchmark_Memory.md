# Benchmarking Screening Methods: Time and Memory

``` r

library(SigBridgeR)
library(Seurat)
```

## Overview

Accuracy and biological signal are only part of the story: a method that
takes hours and tens of gigabytes is often impractical, even when it
performs well.
[`BenchMemory()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemory.md)
runs every configured screening method on the same inputs and records
wall time and memory for each run.

A benchmark is organised in **rounds**. Within one round every method is
screened once. By default the order of methods is shuffled in each
round, so cold caches, lazy package loading and package-level caches are
spread across methods instead of always penalising whichever method
happens to run first.

Each measurement runs in its own background R process created with
[`callr::r_bg()`](https://callr.r-lib.org/reference/r_bg.html). The
inputs are serialised once into a temporary directory and loaded by each
child, so one method cannot leak loaded Python runtimes, package caches
or large result objects into the next measurement.

## Building a configuration

[`BenchMemory()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemory.md)
takes a `BenchMemoryConfig`, which combines the shared benchmark inputs
with the controls returned by
[`BenchMemoryConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfigSet.md):

``` r

tiny_sc <- CreateSeuratObject(
  counts = Matrix::Matrix(
    matrix(
      seq_len(20L),
      nrow = 5L,
      dimnames = list(paste0("Gene", 1:5), paste0("Cell", 1:4))
    ),
    sparse = TRUE
  )
)

config <- BenchMemoryConfig(
  matched_bulk = matrix(
    seq_len(5L),
    nrow = 5L,
    dimnames = list(paste0("Gene", 1:5), "Sample1")
  ),
  sc_data = tiny_sc,
  phenotype = c(Sample1 = 1),
  label_type = "{method_name}",
  phenotype_class = "binary",
  screen_method = "scPP",
  method_param = MethodParamSet(
    scPP_param = list(Log2FC_cutoff = 0.585)
  ),
  bench_config = BenchMemoryConfigSet(
    rounds = 3L,
    order = "random",
    probe_resources = FALSE
  )
)

config@screen_method
#> [1] "scPP"
config@bench_config[c("rounds", "order", "probe_resources")]
#> $rounds
#> [1] 3
#> 
#> $order
#> [1] "random"
#> 
#> $probe_resources
#> [1] FALSE
```

The same structure is used for a real benchmark. The only difference is
the size of the inputs and the number of methods:

``` r

config <- BenchMemoryConfig(
  matched_bulk = bulk,
  sc_data = sc,
  phenotype = phenotype,
  label_type = "{method_name}",
  phenotype_class = "binary",
  screen_method = c("scPP", "scPAS", "LP_SGL"),
  method_param = MethodParamSet(
    scPP_param = list(Log2FC_cutoff = 0.585),
    scPAS_param = list(permutation_times = 500)
  ),
  bench_config = BenchMemoryConfigSet(
    rounds = 3L,
    order = "random",
    measure_rss = TRUE,
    set_threads = 1L,
    timeout_sec = 1800
  )
)
```

`method_param` has one field per registered method and each field is
passed only to the method it is named after. Use
[`MethodParamSet()`](https://wanglabcsu.github.io/sigbridger/reference/MethodParamSet.md)
to get the full set of empty defaults and fill in only the methods you
need.

## Running the benchmark

``` r

res <- BenchMemory(config, verbose = FALSE)
```

The result is a `data.table` with one row per method per round. It
contains the columns returned by
[`bench::mark()`](https://bench.r-lib.org/reference/mark.html) (`min`,
`median`, `itr/sec`, `mem_alloc`, `n_itr`, `n_gc`, `total_time`, …),
plus the benchmark bookkeeping columns:

| Column | Meaning |
|----|----|
| `method_name` | Screening method |
| `round` | Round number |
| `position` | Position of the method within the round |
| `status` | `"ok"`, `"error"` or `"timeout"` |
| `error` | Error message for failed measurements |
| `mem_peak_mb` | Peak R heap high-water mark in MB |
| `rss_peak_mb` | Peak process resident set size in MB (`measure_rss = TRUE`, Linux only) |
| `elapsed` | Wall-clock seconds spent by the measurement |

Successful measurements have `status = "ok"`. Failed or timed-out
measurements are kept with `status = "error"` or `"timeout"`, an `error`
message and `NA` memory metrics, so one broken method does not silently
disappear.

## Summarising the result

Timing is naturally variable, so summarise across rounds rather than
reading a single row. `median` is a `difftime`; convert it to seconds
before aggregating:

``` r

ok <- res[res$status == "ok", ]
ok$median_sec <- as.numeric(ok$median)

aggregate(
  cbind(median_sec, mem_peak_mb, rss_peak_mb) ~ method_name,
  data = ok,
  FUN = median
)
```

The same summary with `data.table`:

``` r

res[
  status == "ok",
  .(
    median_sec = median(as.numeric(median)),
    mem_peak_mb = median(mem_peak_mb),
    rss_peak_mb = median(rss_peak_mb),
    rounds = .N
  ),
  by = method_name
]
```

If a method failed, inspect those rows before drawing conclusions:

``` r

res[status != "ok", c("method_name", "round", "status", "error"), with = FALSE]
```

## Benchmark controls

[`BenchMemoryConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfigSet.md)
collects and validates the controls. The most important ones are:

| Control | Default | Purpose |
|----|----|----|
| `rounds` | `3L` | Number of times every method is screened |
| `order` | `"random"` | Shuffle methods within each round, or keep `"fixed"` order |
| `seed` | package seed | Base seed for shuffling and for each measurement child |
| `timeout_sec` | `10800L` | Maximum seconds for one method in one round |
| `on_failure` | `"skip"` | Record the failure and continue, or `"stop"` the whole benchmark |
| `set_threads` | `1L` | Pin each child process to this many threads |
| `measure_rss` | `FALSE` | Also record peak process RSS (Linux only) |
| `probe_resources` | `TRUE` | Collect a machine resource summary before staging the inputs |
| `min_iterations`, `max_iterations` | `1L`, `1L` | Iteration bounds for [`bench::mark()`](https://bench.r-lib.org/reference/mark.html) |
| `memory`, `filter_gc`, `check`, `min_time`, `relative`, `time_unit` | [`bench::mark()`](https://bench.r-lib.org/reference/mark.html) defaults | Forwarded to [`bench::mark()`](https://bench.r-lib.org/reference/mark.html) |

A few practical choices:

- **Use at least three rounds.** A single round mixes genuine cost with
  one-off startup effects. Odd round counts also let the median ignore
  one unusually slow round.
- **Keep `order = "random"`.** It prevents the first method in the list
  from systematically paying for cold caches. The shuffle is
  reproducible because `seed` plus the round number is used.
- **Pin threads when comparing methods.** `set_threads = 1L` avoids a
  method that internally uses many threads from looking faster merely
  because it consumed more cores. Record the setting alongside the
  results.
- **Set a realistic timeout.** A hung Python backend or an accidental
  full-data copy should not block the benchmark for days.
- **`probe_resources = FALSE`** skips the initial
  [`gc()`](https://rdrr.io/r/base/gc.html) and machine probe, which is
  useful when the workspace is very large and startup time matters.

## Interpreting the numbers

Three memory numbers are reported because they describe different
things:

- `mem_alloc` is the total allocated during the run, i.e. allocation
  churn. A method can allocate a lot and still keep a small peak.
- `mem_peak_mb` is the peak R heap high-water mark. It is reset before
  each measurement and is usually the best indicator of R-side peak
  memory.
- `rss_peak_mb` is the peak resident set size of the whole process. It
  also includes native libraries, Python runtimes and memory mapped
  files, so it is closer to what a scheduler will reserve.

Compare methods measured on the same machine, with the same data, thread
count and package versions. Absolute timings are not portable across
machines, but the relative ranking is often stable.

## See also

- [`BenchMemoryConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfig.md)
  and
  [`BenchMemoryConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfigSet.md)
  for the configuration interface.
- [`BenchAcc()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAcc.md),
  [`BenchGSEA()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEA.md)
  and
  [`BenchGeneScore()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScore.md)
  for accuracy, functional signal and cell-level pathway scores.
- [`bench::mark()`](https://bench.r-lib.org/reference/mark.html) for the
  timing and allocation measurements.
