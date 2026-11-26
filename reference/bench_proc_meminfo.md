# System memory overview

Returns the total and the available physical memory of the machine, in
gigabytes (GB, i.e. 1024^3 bytes).

## Usage

``` r
bench_proc_meminfo()
```

## Value

A named `numeric` vector of length 2 with elements `MemTotal` and
`MemAvailable` (GB). Elements that cannot be determined are `NA_real_`.

## Details

The values are obtained directly from the operating system:

- **Linux**: `MemTotal` and `MemAvailable` from `/proc/meminfo` (falls
  back to `sysconf()` for the total if `/proc` is unreadable;
  `MemAvailable` requires Linux \>= 3.14).

- **macOS**: total from `sysctl(hw.memsize)`; available is approximated
  as free plus inactive pages reported by `host_statistics64()`.

- **Windows**: `ullTotalPhys` and `ullAvailPhys` from
  `GlobalMemoryStatusEx()`.

- **Other Unix-alikes**: total from `sysconf()`; available only if
  `_SC_AVPHYS_PAGES` is supported.

Note that "available" is defined slightly differently by each operating
system, so cross-platform comparisons should be made with care.

## Examples

``` r
if (FALSE) { # \dontrun{
bench_proc_meminfo()
} # }
```
