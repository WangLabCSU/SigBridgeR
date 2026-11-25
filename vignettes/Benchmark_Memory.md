    library(SigBridgeR)
    library(Seurat)

# Overview

Accuracy and biological signal are only part of the story: a method that
takes hours and tens of gigabytes is often impractical, even when it
performs well. `BenchMemory()` runs every configured screening method on
the same inputs and records wall time and memory for each run.

A benchmark is organised in **rounds**. Within one round every method is
screened once. By default the order of methods is shuffled in each
round, so cold caches, lazy package loading and package-level caches are
spread across methods instead of always penalising whichever method
happens to run first.

Each measurement runs in its own background R process created with
`callr::r_bg()`. The inputs are serialised once into a temporary
directory and loaded by each child, so one method cannot leak loaded
Python runtimes, package caches or large result objects into the next
measurement.

# Building a configuration

`BenchMemory()` takes a `BenchMemoryConfig`, which combines the shared
benchmark inputs with the controls returned by `BenchMemoryConfigSet()`:

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
    config@bench_config[c("rounds", "order", "probe_resources")]

The same structure is used for a real benchmark. The only difference is
the size of the inputs and the number of methods:

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

`method_param` has one field per registered method and each field is
passed only to the method it is named after. Use `MethodParamSet()` to
get the full set of empty defaults and fill in only the methods you
need.

# Running the benchmark

    res <- BenchMemory(config, verbose = FALSE)

The result is a `data.table` with one row per method per round. It
contains the columns returned by `bench::mark()` (`min`, `median`,
`itr/sec`, `mem_alloc`, `n_itr`, `n_gc`, `total_time`, …), plus the
benchmark bookkeeping columns:

<table>
<colgroup>
<col style="width: 47%" />
<col style="width: 52%" />
</colgroup>
<thead>
<tr>
<th>Column</th>
<th>Meaning</th>
</tr>
</thead>
<tbody>
<tr>
<td><code>method_name</code></td>
<td>Screening method</td>
</tr>
<tr>
<td><code>round</code></td>
<td>Round number</td>
</tr>
<tr>
<td><code>position</code></td>
<td>Position of the method within the round</td>
</tr>
<tr>
<td><code>status</code></td>
<td><code>"ok"</code>, <code>"error"</code> or
<code>"timeout"</code></td>
</tr>
<tr>
<td><code>error</code></td>
<td>Error message for failed measurements</td>
</tr>
<tr>
<td><code>mem_peak_mb</code></td>
<td>Peak R heap high-water mark in MB</td>
</tr>
<tr>
<td><code>rss_peak_mb</code></td>
<td>Peak process resident set size in MB
(<code>measure_rss = TRUE</code>, Linux only)</td>
</tr>
<tr>
<td><code>elapsed</code></td>
<td>Wall-clock seconds spent by the measurement</td>
</tr>
</tbody>
</table>

Successful measurements have `status = "ok"`. Failed or timed-out
measurements are kept with `status = "error"` or `"timeout"`, an `error`
message and `NA` memory metrics, so one broken method does not silently
disappear.

# Summarising the result

Timing is naturally variable, so summarise across rounds rather than
reading a single row. `median` is a `difftime`; convert it to seconds
before aggregating:

    ok <- res[res$status == "ok", ]
    ok$median_sec <- as.numeric(ok$median)

    aggregate(
      cbind(median_sec, mem_peak_mb, rss_peak_mb) ~ method_name,
      data = ok,
      FUN = median
    )

The same summary with `data.table`:

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

If a method failed, inspect those rows before drawing conclusions:

    res[status != "ok", c("method_name", "round", "status", "error"), with = FALSE]

# Benchmark controls

`BenchMemoryConfigSet()` collects and validates the controls. The most
important ones are:

<table>
<colgroup>
<col style="width: 33%" />
<col style="width: 33%" />
<col style="width: 33%" />
</colgroup>
<thead>
<tr>
<th>Control</th>
<th>Default</th>
<th>Purpose</th>
</tr>
</thead>
<tbody>
<tr>
<td><code>rounds</code></td>
<td><code>3L</code></td>
<td>Number of times every method is screened</td>
</tr>
<tr>
<td><code>order</code></td>
<td><code>"random"</code></td>
<td>Shuffle methods within each round, or keep <code>"fixed"</code>
order</td>
</tr>
<tr>
<td><code>seed</code></td>
<td>package seed</td>
<td>Base seed for shuffling and for each measurement child</td>
</tr>
<tr>
<td><code>timeout_sec</code></td>
<td><code>10800L</code></td>
<td>Maximum seconds for one method in one round</td>
</tr>
<tr>
<td><code>on_failure</code></td>
<td><code>"skip"</code></td>
<td>Record the failure and continue, or <code>"stop"</code> the whole
benchmark</td>
</tr>
<tr>
<td><code>set_threads</code></td>
<td><code>1L</code></td>
<td>Pin each child process to this many threads</td>
</tr>
<tr>
<td><code>measure_rss</code></td>
<td><code>FALSE</code></td>
<td>Also record peak process RSS (Linux only)</td>
</tr>
<tr>
<td><code>probe_resources</code></td>
<td><code>TRUE</code></td>
<td>Collect a machine resource summary before staging the inputs</td>
</tr>
<tr>
<td><code>min_iterations</code>, <code>max_iterations</code></td>
<td><code>1L</code>, <code>1L</code></td>
<td>Iteration bounds for <code>bench::mark()</code></td>
</tr>
<tr>
<td><code>memory</code>, <code>filter_gc</code>, <code>check</code>,
<code>min_time</code>, <code>relative</code>,
<code>time_unit</code></td>
<td><code>bench::mark()</code> defaults</td>
<td>Forwarded to <code>bench::mark()</code></td>
</tr>
</tbody>
</table>

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
- **`probe_resources = FALSE`** skips the initial `gc()` and machine
  probe, which is useful when the workspace is very large and startup
  time matters.

# Interpreting the numbers

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

# See also

- `BenchMemoryConfig()` and `BenchMemoryConfigSet()` for the
  configuration interface.
- `BenchAcc()`, `BenchGSEA()` and `BenchGeneScore()` for accuracy,
  functional signal and cell-level pathway scores.
- `bench::mark()` for the timing and allocation measurements.
