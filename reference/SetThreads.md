# Configure Parallel Execution Backends

Unified interface to configure thread counts and acceleration features
across system-level (OpenMP/data.table/cheapr) and TensorFlow backends.
Uses a hierarchical configuration model where global `threads` serves as
default for all backends unless explicitly overridden.

**Critical**: TensorFlow configuration must be set *before* importing
TensorFlow via `reticulate::import("tensorflow")`.

## Usage

``` r
SetThreads(
  threads = getFuncOption("nthreads"),
  dt = threads,
  cheapr = threads,
  qs2 = threads,
  wgcna = max(threads, 2L),
  seurat = threads,
  openmp = NULL,
  tf_config = list(xla_flag = "--tf_xla_auto_jit=2 --tf_xla_cpu_global_jit", xla_device =
    NULL, inter_op = NULL, intra_op = c(1L, NULL)),
  verbose = getFuncOption("verbose"),
  ...
)
```

## Arguments

- threads:

  Integer. Global thread count used as default for all backends.
  Default: `2L`

- dt:

  Integer. Thread count for data.table (default: inherited from
  `threads`).

- cheapr:

  Integer. Thread count for cheapr (default: inherited from `threads`).

- qs2:

  Integer. Thread count for qs2 (default: inherited from `threads`).

- wgcna:

  Integer. Thread count for WGCNA (default: inherited from
  `max(threads, 2L)`). Values below 2 are raised to 2, because
  [`WGCNA::enableWGCNAThreads()`](https://rdrr.io/pkg/WGCNA/man/allowWGCNAThreads.html)
  requires at least two threads.

- seurat:

  Integer. Thread count for Seurat (default: inherited from `threads`).
  Only applied when Seurat \>= 6.0.0 is installed, as
  [`Seurat::setThreads()`](https://satijalab.org/seurat/reference/setThreads.html)
  is not available in earlier versions.

- openmp:

  Integer. Thread count for OpenMP (default: inherited from `threads`).

- tf_config:

  Named list for TensorFlow-specific configuration:

  - `xla_flag`: Character. XLA JIT compilation flags (default:
    auto-optimized)

  - `xla_device`: Integer. XLA device ID (default: `1L`)

  - `inter_op`: Integer. Inter-op parallelism threads (default:
    `max(2, floor(threads/4))`)

  - `intra_op`: Integer. Intra-op parallelism threads (default: `1L`).
    If `NULL`, inherits global `threads` by default.

- verbose:

  Logical. Whether to print verbose output (default: inherited from
  function options).

- ...:

  Additional arguments for
  [`data.table::setDTthreads()`](https://rdrr.io/pkg/data.table/man/openmp-utils.html)
  (e.g., `restore`).

## Value

Invisible list with old/new values per backend.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage: auto-detect and configure
SetThreads()

# Explicit thread count for CPU-intensive workloads
SetThreads(threads = 12L)

# TensorFlow-optimized configuration for deep learning
SetThreads(
  threads = 8L,
  tf_config = list(
    inter_op = 2L,
    intra_op = 8L
  )
)
# Import AFTER SetThreads()

# Configure only data.table for memory-efficient workflows
SetThreads(dt = 4L , verbose = FALSE)
} # }
```
