# Configure Parallel Execution Backends

Unified interface to configure thread counts and acceleration features
across system-level (OpenMP/data.table) and TensorFlow backends. Uses a
hierarchical configuration model:

- `threads`: Global thread count (default: half physical cores)

- `backend`: System-level backends to configure

- `tf_config`: Named list for TensorFlow-specific settings

## Usage

``` r
setThreads(threads = NULL, backend = c("openmp", "dt"), tf_config = NULL, ...)
```

## Arguments

- threads:

  Integer. Global thread count used as default for all backends. If
  `NULL` (default), uses `floor(availableCores() / 2)`. Applied to:
  OpenMP, data.table, and TensorFlow intra-op (unless overridden).

- backend:

  Character vector. System-level backends to configure: `"openmp"` (sets
  `OMP_NUM_THREADS`), `"dt"` (data.table threads). Default:
  `c("openmp", "dt")`.

- tf_config:

  Named list for TensorFlow-specific configuration:

  - `xla`: Logical. Enable XLA JIT compilation (default: `FALSE`)

  - `inter_op`: Integer. Inter-op parallelism threads (default:
    auto-derived)

  - `intra_op`: Integer. Intra-op parallelism threads (default:
    `threads`)

  **Must be set BEFORE importing TensorFlow**.

- ...:

  Additional arguments for
  [`data.table::setDTthreads()`](https://rdrr.io/pkg/data.table/man/openmp-utils.html)
  (e.g., `restore`). Includes `verbose` logical flag (default: `TRUE`).

## Value

Invisible list with old/new values per backend.

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal: auto-configure all backends with half cores
setThreads()

# Explicit global thread count
setThreads(threads = 8)

# TensorFlow optimization (MUST run BEFORE importing TensorFlow)
setThreads(
  threads = 8,
  tf_config = list(
    xla = TRUE,
    inter_op = 2,
    intra_op = 8
  )
)
reticulate::import("tensorflow")  # Import AFTER configuration

# Fine-grained control: override specific backends
setThreads(
  threads = 16,
  backend = "dt",  # Only configure data.table
  tf_config = list(inter_op = 4)  # intra_op inherits threads=16
)
} # }
```
