# Resolve Benchmark Controls

Merges the user supplied `bench_config` with the defaults and validates
the names, so that a typo aborts instead of being silently ignored.

## Usage

``` r
BenchMemoryConfigSet(
  rounds = 3L,
  probe_resources = TRUE,
  min_iterations = 1L,
  max_iterations = 1L,
  order = c("random", "fixed"),
  seed = getFuncOption("seed"),
  memory = TRUE,
  filter_gc = FALSE,
  check = FALSE,
  measure_rss = FALSE,
  on_failure = c("skip", "stop"),
  set_threads = 1L,
  min_time = 0.5,
  relative = FALSE,
  time_unit = NULL,
  timeout_sec = 10800L,
  ...
)
```

## Arguments

- rounds:

  `integer`. Number of benchmark rounds; every method is screened once
  per round and per-method summaries are taken across rounds. Odd values
  let the median ignore a single cold round. Must be at least 1.
  Default: `3L`.

- probe_resources:

  `logical`. Whether
  [`BenchMemory()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemory.md)
  should call
  [`bench_probe_resources()`](https://wanglabcsu.github.io/sigbridger/reference/bench_probe_resources.md)
  before staging the inputs. The probe calls
  [`gc()`](https://rdrr.io/r/base/gc.html) and can be expensive on large
  workspaces; disable it when the machine resource summary is not
  needed. Default: `TRUE`.

- min_iterations:

  Each expression will be evaluated a minimum of `min_iterations` times.

- max_iterations:

  Each expression will be evaluated a maximum of `max_iterations` times.

- order:

  `character`. Order in which methods run within a round: `"random"`
  shuffles them (seeded by `seed`) to spread cold-cache and
  package-caching effects; `"fixed"` keeps the given method order.
  Default: `c("random", "fixed")` (first value used).

- seed:

  `integer`. Base seed for shuffling method order each round; the round
  number is added to it. Also set inside each measurement child process.
  Default: `123L`.

- memory:

  If `TRUE` (the default when R is compiled with memory profiling),
  track memory allocations using
  [`utils::Rprofmem()`](https://rdrr.io/r/utils/Rprofmem.html). If
  `FALSE` disable memory tracking.

- filter_gc:

  If `TRUE` remove iterations that contained at least one garbage
  collection before summarizing. If `TRUE` but an expression had a
  garbage collection in every iteration, filtering is disabled, with a
  warning.

- check:

  Check if results are consistent. If `TRUE`, checking is done with
  [`all.equal()`](https://rdrr.io/r/base/all.equal.html), if `FALSE`
  checking is disabled and results are not stored. If `check` is a
  function that function will be called with each pair of results to
  determine consistency.

- measure_rss:

  `logical`. Whether to record peak process resident set size (RSS) in
  addition to R heap memory. Supported on Linux only. Default: `FALSE`.

- on_failure:

  `character`. Action when a method times out or errors in a round:
  `"skip"` warns and continues with the remaining methods, `"stop"`
  aborts the whole benchmark. Default: `c("skip", "stop")` (first value
  used).

- set_threads:

  `integer`. Number of threads the measurement child process is pinned
  to via
  [`SetThreads()`](https://wanglabcsu.github.io/sigbridger/reference/SetThreads.md),
  so runs stay comparable across methods. Default: `1L`.

- min_time:

  The minimum number of seconds to run each expression, set to `Inf` to
  always run `max_iterations` times instead.

- relative:

  If `TRUE` all summaries are computed relative to the minimum execution
  time rather than absolute time.

- time_unit:

  If `NULL` the times are reported in a human readable fashion depending
  on each value. If one of 'ns', 'us', 'ms', 's', 'm', 'h', 'd', 'w' the
  time units are instead expressed as nanoseconds, microseconds,
  milliseconds, seconds, hours, minutes, days or weeks respectively.

- timeout_sec:

  `numeric`. Maximum wall-clock seconds allowed for one method in one
  round; the child process is killed when the limit is exceeded. Use
  `Inf` to disable the limit. Default: `10800L` (3 hours).

- ...:

  Expressions to benchmark, if named the `expression` column will be the
  name, otherwise it will be the deparsed expression.

## Value

A named `list` of resolved controls.

## See also

Other ConfigSet:
[`BenchGSEAConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEAConfigSet.md),
[`BenchGeneScoreConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfigSet.md),
[`MethodParamSet()`](https://wanglabcsu.github.io/sigbridger/reference/MethodParamSet.md)

## Examples

``` r
# \donttest{
BenchMemoryConfigSet()
#> $rounds
#> [1] 3
#> 
#> $probe_resources
#> [1] TRUE
#> 
#> $min_iterations
#> [1] 1
#> 
#> $max_iterations
#> [1] 1
#> 
#> $order
#> [1] "random"
#> 
#> $seed
#> [1] 123
#> 
#> $memory
#> [1] TRUE
#> 
#> $filter_gc
#> [1] FALSE
#> 
#> $check
#> [1] FALSE
#> 
#> $measure_rss
#> [1] FALSE
#> 
#> $on_failure
#> [1] "skip"
#> 
#> $set_threads
#> [1] 1
#> 
#> $min_time
#> [1] 0.5
#> 
#> $relative
#> [1] FALSE
#> 
#> $time_unit
#> NULL
#> 
#> $timeout_sec
#> [1] 10800
#> 
# }
```
