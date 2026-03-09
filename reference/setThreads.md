# Configure Parallel Execution Backends

Unified interface to configure thread counts and acceleration features
across system-level (OpenMP/data.table/cheapr) and TensorFlow backends.
Uses a hierarchical configuration model where global `threads` serves as
default for all backends unless explicitly overridden.

**Critical**: TensorFlow configuration must be set *before* importing
TensorFlow via `reticulate::import("tensorflow")`.

## Usage

``` r
setThreads(
  threads = NULL,
  backend = c("openmp", "dt", "cheapr"),
  tf_config = list(xla_flag = "--tf_xla_auto_jit=2 --tf_xla_cpu_global_jit", xla_device =
    NULL, inter_op = NULL, intra_op = c(1L, NULL)),
  verbose = getFuncOption("verbose"),
  ...
)
```

## Arguments

- threads:

  Integer. Global thread count used as default for all backends. If
  `NULL` (default), uses `floor(availableCores() / 2)`. Applied to:
  OpenMP, data.table, and TensorFlow intra-op (unless overridden).

- backend:

  Character vector. System-level backends to configure: `"openmp"` (sets
  `OMP_NUM_THREADS`), `"dt"` (data.table threads). Default:
  `c("openmp", "dt", "cheapr")`.

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
setThreads()

# Explicit thread count for CPU-intensive workloads
setThreads(threads = 12L)

# TensorFlow-optimized configuration for deep learning
setThreads(
  threads = 8L,
  tf_config = list(
    inter_op = 2L,
    intra_op = 8L
  )
)
library(tensorflow)  # Import AFTER setThreads()

# Configure only data.table for memory-efficient workflows
setThreads(threads = 4, backend = "dt", verbose = FALSE)
} # }
```
