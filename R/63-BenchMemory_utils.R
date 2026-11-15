#' @title Probe Machine Resources
#' @description
#' Collects CPU count, configured thread count, memory figures and load
#' average. Values that cannot be detected are returned as `NA` instead of
#' being guessed.
#' @return A named `list`.
#' @keywords internal
bench_probe_resources <- function() {
  rlang::check_installed("parallel")
  mem <- bench_proc_meminfo() # cpp

  list(
    cores = parallel::detectCores(),
    nthreads = getFuncOption("nthreads"),
    mem_total_mb = mem[["MemTotal"]],
    mem_avail_mb = mem[["MemAvailable"]],
    load1 = bench_proc_load(), # cpp
    heap_used_mb = bench_heap_used_mb()
  )
}


bench_heap_used_mb <- function() {
  sum(gc(verbose = FALSE)[, 6L])
}


bench_tmp_dir <- function(sc_data, bulk_data, pheno_data) {
  tmp_dir <- file.path(tempfile(), paste0("SigBridgeR_bench_tmp"))
  dir.create(tmp_dir, recursive = TRUE)
  n_cores <- max(floor(parallel::detectCores() / 2L), 1L)
  qs2::qs_save(sc_data, file.path(tmp_dir, "sc_data.qs2"), nthreads = n_cores)
  qs2::qd_save(bulk_data, file.path(tmp_dir, "matched_bulk.qdata"), nthreads = n_cores)
  qs2::qd_save(pheno_data, file.path(tmp_dir, "phenotype.qdata"), nthreads = n_cores)
  tmp_dir
}


#' @title Collect Child Process Output
#' @description
#' Appends the buffered `stdout`/`stderr` of a [callr::r_bg()] child to a rolling
#' log. Only the last `n` lines are kept, so a chatty child cannot grow the
#' parent heap during a long measurement.
#' @param log
#'
#' @return A `character` vector of at most `n` lines.
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
#' Runs `method_name` once, optionally after an untimed warm-up, inside
#' [bench::mark()] and returns a single-row table with the timing and memory
#' metrics plus the elapsed wall time of the whole measurement.
#' @keywords internal
bench_mem_measure <- function(config, method_name, round, position, tmp_dir, verbose = FALSE) {
  gc(reset = TRUE)
  started <- Sys.time()

  method_config <- config@method_param[[method_name]]

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
      SigBridgeR::setThreads(bench_config$set_threads, verbose = FALSE)

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
      rss_peak <- if (isTRUE(config@bench_config$measure_rss)) bench_peak_rss_mb() else NA_real_
      list(bench_res = bench_res, mem_peak = mem_peak, rss_peak = rss_peak)
    },
    args = list(
      tmp_dir = tmp_dir,
      n_cores = max(floor(parallel::detectCores() / 2L), 1L),
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
      cli::cli_warn("Cleaning up the child process for {method_name} in round {round}")
      if (bench_process$is_alive()) {
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
      error = error_msg,
      started = started
    ))
  }

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
        tips = "The benchmark child process died; see the returned {.code error} column"
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

  gc(verbose = FALSE)

  res
}

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
