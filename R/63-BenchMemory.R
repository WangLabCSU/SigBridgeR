#' @title Benchmark Screening Methods: Wall Time and Memory
#'
#' @description
#' Runs every screening method listed in a [BenchMemoryConfig] on the same data
#' and measures its wall time and memory footprint with [bench::mark()].
#'
#' @details
#' The benchmark is organised in *rounds*. Within a round every method is
#' screened once, and the round order is shuffled by default so that cold-cache,
#' lazy-loading and package-level caching effects are spread across methods
#' instead of always penalising the first one. One row is returned per method
#' per round; summaries across rounds are left to the caller.
#'
#' Each measurement runs in its own background R process ([callr::r_bg()]) with
#' its own copies of the inputs, staged in a temporary directory with `qs2`.
#' A child process therefore cannot leak package-level caches, loaded Python
#' runtimes or large result objects into the next measurement. The child is
#' always terminated on exit so it cannot outlive the call, and `timeout_sec`
#' bounds the wall time of a single method in a single round.
#'
#' Three memory numbers are reported because they mean different things:
#'
#' * `mem_alloc` - bytes allocated during the run (allocation churn, from
#'   [bench::mark()]).
#' * `mem_peak_mb` - peak R heap high-water mark, from `gc(reset = TRUE)`.
#' * `rss_peak_mb` - peak process resident set size, recorded only when
#'   `measure_rss = TRUE` (Linux only).
#'
#' @section Benchmark controls (`bench_config`):
#' * `rounds` (`3L`): number of rounds; every method is screened once per
#'   round.
#' * `min_iterations`, `max_iterations` (`1L`, `1L`): iteration bounds passed
#'   to [bench::mark()].
#' * `order` (`"random"`): `"random"` shuffles methods within each round using
#'   `seed` plus the round number; `"fixed"` keeps the given order.
#' * `seed` (`getFuncOption("seed")`): base seed for shuffling and for each
#'   measurement child.
#' * `memory`, `filter_gc`, `check`, `min_time`, `relative`, `time_unit`:
#'   forwarded to [bench::mark()].
#' * `measure_rss` (`FALSE`): also record peak process RSS (Linux only).
#' * `on_failure` (`"skip"`): `"skip"` records the failure and continues,
#'   `"stop"` aborts the whole benchmark.
#' * `set_threads` (`1L`): number of threads the child process is pinned to.
#' * `timeout_sec` (`10800L`): maximum wall-clock seconds for one method in one
#'   round; use `Inf` to disable the limit.
#'
#' @param obj An object of class [BenchMemoryConfig]
#' @param ... Additional arguments passed to methods. Currently unused.
#' @param verbose `logical`. Whether to report progress.
#'   Default: `getFuncOption("verbose")`.
#'
#' @return A `data.table` with one row per screening method per round. Each row
#'   contains the [bench::mark()] columns (`min`, `median`, `itr/sec`,
#'   `mem_alloc`, `n_itr`, `n_gc`, `total_time`, ...), with `expression`
#'   replaced by the method name, plus `method_name`, `round`, `position`,
#'   `status`, `error`, `mem_peak_mb`, `rss_peak_mb` and `elapsed`. Successful
#'   measurements have `status = "ok"`; failed or timed-out measurements are
#'   kept with `status` set to `"error"` or `"timeout"` and `NA` memory
#'   metrics, while their benchmark columns are filled with `NA`.
#'
#' @examples
#' \dontrun{
#' config <- BenchMemoryConfig(
#'   matched_bulk = bulk,
#'   sc_data = sc,
#'   phenotype = pheno,
#'   screen_method = c("scPP", "scPAS", "LP_SGL"),
#'   method_param = MethodParamSet(list(
#'     scPP_param = list(Log2FC_cutoff = 0.585),
#'     scPAS_param = list(permutation_times = 500)
#'   )),
#'   bench_config = list(rounds = 3L, timeout_sec = 1800)
#' )
#' res <- BenchMemory(config)
#' res
#' res[status == "ok", .(median_sec = median(median)), by = method_name]
#' }
#'
#' @seealso [BenchMemoryConfig()], [BenchMemoryConfigSet()], [bench::mark()]
#' @family Benchmark
#'
#' @name BenchMemory
#' @export
BenchMemory <- function(obj, ...) {
  UseMethod("BenchMemory")
}

#' @rdname BenchMemory
#' @export
`BenchMemory.SigBridgeR::BenchMemoryConfig` <- function(
  obj,
  ...,
  verbose = getFuncOption("verbose")
) {
  rlang::check_installed(c("bench", "glue", "callr", "qs2"))

  methods <- unique(obj@screen_method)

  if (verbose) {
    cli::cli_progress_step("Detecting machine resources")
  }

  resources <- bench_probe_resources()

  if (verbose) {
    cli::cli_progress_step("Initializing temporary directory")
  }

  # * seed & threads & init
  tmp_dir <- bench_tmp_dir(obj@sc_data, obj@matched_bulk, obj@phenotype)
  on.exit(
    {
      cli::cli_progress_done()
      unlink(tmp_dir, recursive = TRUE, force = TRUE)
      gc(verbose = FALSE)
    },
    add = TRUE
  )

  rows <- vector("list", length(methods) * obj@bench_config$rounds)
  rounds <- obj@bench_config$rounds

  for (round in seq_len(rounds)) {
    order <- if (obj@bench_config$order == "random") {
      set.seed(obj@bench_config$seed + round)
      sample(methods)
    } else {
      methods
    }

    for (i in seq_along(order)) {
      method_name <- order[[i]]

      if (isTRUE(verbose)) {
        cli::cli_progress_step(
          "Round {round}/{rounds}: {.val {method_name}}"
        )
      }

      rows[[(round - 1L) * length(methods) + i]] <- bench_mem_measure(
        config = obj,
        method_name = method_name,
        round = round,
        position = i,
        tmp_dir = tmp_dir,
        verbose = verbose
      )
    }
  }

  metrics <- data.table::rbindlist(rows, fill = TRUE)
  data.table::setorder(metrics, c("method_name", "round", "position"))
  cli::cli_progress_done()

  metrics
}

#' @rdname BenchMemory
#' @export
BenchMemory.default <- function(obj, ...) {
  cls_obj <- class(obj)
  avaialble_methods <- utils::methods(BenchMemory, all.names = TRUE) |>
    as.character() |>
    gsub(pattern = ".*\\.", replacement = "", x = _)
  Abort(
    "BenchMemory is not implemented for class: {.cls {cls_obj}}",
    "Available methods: {.cls {avaialble_methods}}"
  )
}
