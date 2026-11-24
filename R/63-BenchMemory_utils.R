#' @title Probe Machine Resources
#' @description
#' Collects CPU count, configured thread count, physical memory, load average
#' and R heap usage. Values that cannot be detected are returned as `NA`
#' instead of being guessed.
#'
#' @return A named `list` with elements:
#'   * `cores`: CPU count from `parallel::detectCores()` (may be `NA`);
#'   * `nthreads`: configured thread count from `getFuncOption("nthreads")`;
#'   * `mem_total_gb`, `mem_avail_gb`: total and available physical memory as
#'     reported by [bench_proc_meminfo()];
#'   * `load1`: 1-minute load average from [bench_proc_load()] (`NA_real_` on
#'     platforms without a load average);
#'   * `heap_used_mb`: R heap high-water mark from `bench_heap_used_mb()`.
#'
#' @keywords internal
bench_probe_resources <- function() {
  rlang::check_installed("parallel")
  mem <- bench_proc_meminfo() # cpp

  list(
    cores = parallel::detectCores(),
    nthreads = getFuncOption("nthreads"),
    mem_total_gb = mem[["MemTotal"]],
    mem_avail_gb = mem[["MemAvailable"]],
    load1 = bench_proc_load(), # cpp
    heap_used_mb = bench_heap_used_mb()
  )
}


#' @title R Heap High-Water Mark
#' @description
#' Returns the peak size of the R heap in megabytes, taken from the `max used`
#' figures of `gc()`. A `gc(reset = TRUE)` before a measurement resets this
#' high-water mark, so later calls report the peak reached during that
#' measurement rather than since the process started.
#'
#' @return A single `numeric` value: peak R heap usage in MB.
#'
#' @keywords internal
bench_heap_used_mb <- function() {
  sum(gc(verbose = FALSE)[, 6L])
}


#' @title Stage Benchmark Inputs in a Temporary Directory
#' @description
#' Serialises the three shared benchmark inputs into a fresh temporary
#' directory with `qs2`, so that each measurement child process can load its
#' own copies instead of sharing mutable objects with the parent.
#'
#' @param sc_data A `Seurat` object holding the cells to screen.
#' @param bulk_data The matched bulk data passed to [Screen()].
#' @param pheno_data The phenotype data passed to [Screen()].
#'
#' @details Serialisation uses half of the detected CPU cores (at least one)
#'   as the `qs2` thread count.
#'
#' @return The path of the created directory, containing `sc_data.qs2`,
#'   `matched_bulk.qdata` and `phenotype.qdata` and the number of CPU cores.
#'   The caller is responsible for removing it, as [BenchMemory()] does via `on.exit()`.
#'
#' @keywords internal
bench_tmp_dir <- function(sc_data, bulk_data, pheno_data) {
  tmp_dir <- file.path(tempfile(), paste0("SigBridgeR_bench_tmp"))
  dir.create(tmp_dir, recursive = TRUE)
  n_cores <- max(floor(parallel::detectCores() / 2L), 2L) %||% 2L

  qs2::qs_save(sc_data, file.path(tmp_dir, "sc_data.qs2"), nthreads = n_cores, compress_level = 1L)
  qs2::qd_save(bulk_data, file.path(tmp_dir, "matched_bulk.qdata"), nthreads = n_cores, compress_level = 1L)
  qs2::qd_save(pheno_data, file.path(tmp_dir, "phenotype.qdata"), nthreads = n_cores, compress_level = 1L)

  list(tmp_dir = tmp_dir, n_cores = n_cores)
}


#' @title Collect Child Process Output
#' @description
#' Appends the buffered `stdout`/`stderr` of a [callr::r_bg()] child to a rolling
#' log. Only the last `n` lines are kept, so a chatty child cannot grow the
#' parent heap during a long measurement.
#' @param log `character`. Lines collected by earlier calls; the returned
#'   vector keeps at most `n` of the combined old and new lines.
#' @param process A [callr::r_bg()] handle whose output and error pipes are
#'   read. Pipe errors (e.g. after the child has been reaped) are ignored, so
#'   already collected lines are preserved.
#' @param n `integer`. Maximum number of lines to keep. Default: `100L`.
#'
#' @return A `character` vector of at most `n` lines.
#'
#' @keywords internal
bench_child_log <- function(log, process, n = 100L) {
  # * processx closes the pipes once the child is reaped, so a late read raises
  # * instead of returning nothing: keep whatever was already collected
  new <- tryCatch(
    c(process$read_output(), process$read_error()),
    error = function(e) character()
  )

  utils::tail(c(log, new), n)
}

#' @title Measure One Screening Method
#' @description
#' Runs `method_name` once inside [bench::mark()] in an isolated
#' [callr::r_bg()] child and returns a one-row table with the timing and
#' memory metrics. The parent supervises the child: its pipes are drained
#' while it runs, and it is killed when `bench_config$timeout_sec` is exceeded.
#'
#' @param config A [BenchMemoryConfig] object supplying `label_type`,
#'   `phenotype_class` and the `bench_config` controls (`seed`, `set_threads`,
#'   `measure_rss`, `on_failure`, `timeout_sec`, and the arguments forwarded to
#'   [bench::mark()]).
#' @param method_name `character`. Name of one screening method; its arguments
#'   are taken from `config@method_param[[paste0(method_name, "_param")]]`.
#' @param round `integer`. Round number, recorded in the result.
#' @param position `integer`. Position of the method within the round.
#' @param tmp_dir `character`. Directory created by `bench_tmp_dir()`` holding
#'   the serialised inputs loaded by the child.
#' @param verbose `logical`. Whether to report progress. Default: `FALSE`.
#' @param n_cores `integer`. Number of CPU cores to configure for the child. Used not in
#'   the child, but in the parent.
#'
#' @details
#' The child sets the configured seed and threads, loads the inputs from
#' `tmp_dir`, and runs [Screen()] through [bench::mark()]. It returns the mark
#' result together with the R heap peak (`bench_heap_used_mb()`) and, when
#' `measure_rss` is enabled, the process RSS peak (`bench_peak_rss_mb()``).
#'
#' Failure and timeout are handled according to `on_failure`: `"stop"` aborts,
#' while `"skip"` returns a row from [fail_row()] with `status` set to
#' `"error"` or `"timeout"`. The child is always killed and reaped on exit so
#' it cannot outlive the call.
#'
#' @return A one-row `data.table`: the [bench::mark()] columns, with
#'   `expression` replaced by `method_name`, plus `method_name`, `round`,
#'   `position`, `status`, `error`, `mem_peak_mb`, `rss_peak_mb` and `elapsed`
#'   (total wall-clock seconds of the child run).
#'
#' @keywords internal
bench_mem_measure <- function(
  config,
  method_name,
  round,
  position,
  tmp_dir,
  verbose = FALSE,
  n_cores = 2L
) {
  started <- Sys.time()

  timeout_sec <- config@bench_config$timeout_sec
  # * absolute deadline, `NULL` means "no limit" (`timeout_sec = Inf`)
  deadline <- if (is.finite(timeout_sec)) started + timeout_sec else Inf

  bench_process <- callr::r_bg(
    func = function(
      tmp_dir,
      n_cores,
      label_type,
      phenotype_class,
      bench_config,
      method_config,
      method_name,
      verbose,
      bench_peak_rss_mb,
      bench_heap_used_mb
    ) {
      # * seed & threads
      set.seed(bench_config$seed)
      SigBridgeR::SetThreads(bench_config$set_threads, verbose = FALSE)

      sc_data <- qs2::qs_read(file.path(tmp_dir, "sc_data.qs2"), nthreads = n_cores)
      matched_bulk <- qs2::qd_read(file.path(tmp_dir, "matched_bulk.qdata"), nthreads = n_cores)
      phenotype <- qs2::qd_read(file.path(tmp_dir, "phenotype.qdata"), nthreads = n_cores)

      bench_res <- tryCatch(
        rlang::exec(
          .fn = bench::mark,
          run = {
            res_list <- SigBridgeR::Screen(
              matched_bulk = matched_bulk,
              sc_data = sc_data,
              phenotype = phenotype,
              label_type = label_type,
              phenotype_class = phenotype_class,
              screen_method = method_name,
              !!!method_config
            )
          },
          !!!SigBridgeRUtils::FilterArgs4Func(args_list = bench_config, fun = bench::mark)
        ),
        error = function(e) e
      )
      mem_peak <- bench_heap_used_mb()
      rss_peak <- if (isTRUE(bench_config$measure_rss)) bench_peak_rss_mb() else NA_real_
      list(bench_res = bench_res, mem_peak = mem_peak, rss_peak = rss_peak)
    },
    args = list(
      tmp_dir = tmp_dir,
      n_cores = n_cores,
      label_type = glue::glue(config@label_type), # `method_name` exists
      phenotype_class = config@phenotype_class,
      bench_config = config@bench_config,
      method_config = config@method_param[[paste0(method_name, "_param")]],
      method_name = method_name,
      verbose = verbose,
      bench_peak_rss_mb = bench_peak_rss_mb,
      bench_heap_used_mb = bench_heap_used_mb
    )
  )
  # * the child owns a whole R session (and, for some methods, a python runtime):
  # * it must never outlive this call, not even when the benchmark is interrupted
  # * from the console
  on.exit(
    {
      if (bench_process$is_alive()) {
        cli::cli_warn("Cleaning up the child process for {method_name} in round {round}")
        invisible(bench_process$kill_tree())
      }
    },
    add = TRUE
  )

  # * supervise instead of blocking in `get_result()`, which has no timeout. The
  # * pipes are drained on every iteration: a chatty child (cli, Seurat, python)
  # * would otherwise block on a full pipe buffer, which is indistinguishable
  # * from a method that hangs.
  child_log <- character()
  timed_out <- FALSE

  repeat {
    # update log
    child_log <- bench_child_log(log = child_log, process = bench_process)

    if (!bench_process$is_alive()) {
      break
    }

    if (Sys.time() >= deadline) {
      cli::cli_warn(
        "Child process method {method_name} in round {round} exceeded the timeout({timeout_sec}s) and was killed"
      )
      timed_out <- TRUE
      break
    }

    bench_process$poll_io(250L)
  }

  # --------------------------------------------------------------------------------------------
  # Timeout
  # * collapse the child output into one line: it goes into `error` and into cli
  # * messages, where newlines would break the layout
  child_log <- bench_child_log(log = child_log, process = bench_process)
  child_tail <- gsub("[[:space:]]+", " ", trimws(paste(child_log, collapse = " ")))
  child_tail <- substr(child_tail, 1L, 300L)

  log_suffix <- if (nzchar(child_tail)) {
    paste0("; last child output: ", child_tail)
  } else {
    "; the child printed nothing"
  }

  if (timed_out) {
    cli::cli_warn("Cleaning up the child process for {method_name} in round {round}")

    invisible(bench_process$kill_tree())
    bench_process$wait(5000L) # reap the killed child

    if (identical(config@bench_config$on_failure, "stop")) {
      Abort(
        "Method {.val {method_name}} timed out in round {round}",
        tips = "Check the input data, raise {.code bench_config$timeout_sec}, or set it to {.val Inf} to disable the limit"
      )
    }

    if (nzchar(child_tail)) {
      cli::cli_alert_info("Last child output: {child_tail}")
    }

    return(fail_row(
      method_name = method_name,
      round = round,
      position = position,
      status = "timeout",
      error = "timeout",
      started = started
    ))
  }
  # ----------------------------------------------------------------------------------------------
  # Failure
  # * a child that dies before returning (segfault, OOM kill, `q()`) raises here
  # * instead of returning a value; treat it as a regular failure so that
  # * `on_failure = "skip"` still keeps the remaining methods
  child <- tryCatch(bench_process$get_result(), error = function(e) e)

  if (inherits(child, "error")) {
    error_msg <- sprintf(
      "child process died before returning a result: %s%s",
      conditionMessage(child),
      log_suffix
    )

    if (identical(config@bench_config$on_failure, "stop")) {
      Abort(
        "Method {.val {method_name}} failed in round {round}",
        "The benchmark child process died; see the returned {.code error} column"
      )
    }

    cli::cli_alert_danger("{.val {method_name}} failed in round {round}")
    if (nzchar(child_tail)) {
      cli::cli_text("Last child output: {child_tail}")
    }

    return(fail_row(
      method_name = method_name,
      round = round,
      position = position,
      status = "error",
      error = error_msg,
      started = started
    ))
  }

  # -----------------------------------------------------------------------------------------------
  # success

  bench_res <- mem_peak <- rss_peak <- NULL # ease check NOTE
  c(bench_res, mem_peak, rss_peak) %<-% child

  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  if (inherits(bench_res, "error")) {
    if (identical(config@bench_config$on_failure, "stop")) {
      Abort(
        "Method {.val {method_name}} failed in round {round}",
        tips = conditionMessage(bench_res)
      )
    }
    cli::cli_alert_danger(c(
      "{.val {method_name}} failed in round {round}:",
      "{conditionMessage(bench_res)}"
    ))

    return(fail_row("error", conditionMessage(bench_res)))
  }

  res <- data.table::as.data.table(bench_res)
  if ("expression" %chin% names(res)) {
    data.table::set(res, j = "expression", value = method_name)
  }

  res[, `:=`(
    method_name = method_name,
    round = round,
    position = position,
    status = "ok",
    error = NA_character_,
    mem_peak_mb = mem_peak,
    rss_peak_mb = rss_peak,
    elapsed = elapsed
  )]

  res
}

#' @title Build a Failed-Measurement Row
#' @description
#' Creates the one-row table returned by [bench_mem_measure()] when a child
#' times out or fails, keeping the same columns as a successful row but with
#' `NA` memory metrics.
#'
#' @param method_name `character`. Name of the failed screening method.
#' @param round,position `integer`. Round number and position within the round.
#' @param status `character`. Failure status, e.g. `"error"` or `"timeout"`.
#' @param error `character`. Message describing the failure.
#' @param started `POSIXct` time at which the measurement started, used to
#'   compute `elapsed`.
#'
#' @return A one-row `data.table`.
#'
#' @keywords internal
fail_row <- function(method_name, round, position, status, error, started) {
  data.table::data.table(
    method_name = method_name,
    round = round,
    position = position,
    status = status,
    error = error,
    mem_peak_mb = NA_real_,
    rss_peak_mb = NA_real_,
    elapsed = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )
}
