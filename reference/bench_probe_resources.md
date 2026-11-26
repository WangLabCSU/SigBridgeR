# Probe Machine Resources

Collects CPU count, configured thread count, physical memory, load
average and R heap usage. Values that cannot be detected are returned as
`NA` instead of being guessed.

## Usage

``` r
bench_probe_resources()
```

## Value

A named `list` with elements:

- `cores`: CPU count from
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
  (may be `NA`);

- `nthreads`: configured thread count from `getFuncOption("nthreads")`;

- `mem_total_gb`, `mem_avail_gb`: total and available physical memory as
  reported by
  [`bench_proc_meminfo()`](https://wanglabcsu.github.io/sigbridger/reference/bench_proc_meminfo.md);

- `load1`: 1-minute load average from
  [`bench_proc_load()`](https://wanglabcsu.github.io/sigbridger/reference/bench_proc_load.md)
  (`NA_real_` on platforms without a load average);

- `heap_used_mb`: R heap high-water mark from
  [`bench_heap_used_mb()`](https://wanglabcsu.github.io/sigbridger/reference/bench_heap_used_mb.md).
