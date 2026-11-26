# Peak resident memory of the current R process

Returns the peak resident set size (RSS) reached by the current R
process since it started, in megabytes (MB, i.e. 1024^2 bytes).

## Usage

``` r
bench_peak_rss_mb()
```

## Value

A single `numeric` value: the peak RSS in MB, or `NA_real_` if it cannot
be determined on the current platform.

## Details

The value is obtained directly from the operating system:

- **Windows**: `PeakWorkingSetSize` from `GetProcessMemoryInfo()`.

- **macOS**: `ru_maxrss` from `getrusage(RUSAGE_SELF)`, reported in
  bytes by the kernel.

- **Linux and other Unix-alikes**: `ru_maxrss` from
  `getrusage(RUSAGE_SELF)`, reported in kilobytes by the kernel
  (equivalent to `VmHWM` in `/proc/self/status`).

The peak is monotonically non-decreasing over the lifetime of the
process; it is not reset by garbage collection.

## Examples

``` r
if (FALSE) { # \dontrun{
bench_peak_rss_mb()

# Allocate roughly 400 MB and observe the peak grow
x <- numeric(5e7)
bench_peak_rss_mb()
} # }
```
