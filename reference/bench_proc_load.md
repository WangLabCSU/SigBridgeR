# System load average

Returns the 1-minute system load average.

## Usage

``` r
bench_proc_load()
```

## Value

A single `numeric` value: the 1-minute load average, or `NA_real_` if it
is unavailable on the current platform.

## Details

On Linux, macOS and other Unix-alikes the value comes from
`getloadavg()` (on Linux this is the first field of `/proc/loadavg`).
Windows has no load-average concept, so `NA_real_` is returned there.
Inside containers the value usually reflects the host, not the
container.

## Examples

``` r
if (FALSE) { # \dontrun{
bench_proc_load()
} # }
```
