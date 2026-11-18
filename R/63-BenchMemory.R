#' @title Benchmark Screening Methods: Wall Time and Memory
#'
#' @description
#' Runs every screening method listed in a [BenchMemoryConfig] on the same data
#' and measures its wall time and memory footprint with [bench::mark()].
#'
#' @details
#' The benchmark is organised in *rounds*. Within a round every method is
#' screened once, and the round order is shuffled so that cold-cache,
#' lazy-loading and package-level caching effects are spread across methods
#' instead of always penalising the first one. Per-method summaries are taken
#' across rounds (median, min, max, MAD, CV) rather than from a single run.
#'
#' Measurement is deliberately one `bench::mark()` call per method per round: a
#' single call over all methods would keep every `ScreenMethodResult` alive in
#' the same result object, inflating the heap seen by the methods measured
#' afterwards. Each call drops its result immediately and is followed by `gc()`.
#'
#' Two memory numbers are reported because they mean different things:
#'
#' * `mem_alloc_mb` - bytes allocated during the run (allocation churn, from
#'   [bench::mark()]).
#' * `mem_peak_mb` - peak R heap high-water mark, from `gc(reset = TRUE)`.
#'
#' @section Machine protection:
#' By default (`auto_budget = TRUE`) the benchmark inspects the machine and
#' adapts itself:
#'
#' * threads are pinned to 1 while measuring and restored afterwards, because
#'   multi-threaded runs are neither comparable across methods nor
#'   representative of an interactive session;
#' * with more than one round, the elapsed time of finished rounds is used to
#'   project the total cost, and the remaining rounds are dropped when
#'   `time_budget_sec` would be exceeded;
#' * before each round the free memory is compared with `mem_budget_mb`; a run
#'   that cannot fit is stopped instead of pushing the machine into swap.
#'
#' Every decision is recorded in `attr(result, "settings")`, and all thresholds
#' can be overridden through `bench_config` of the configuration.
#'
#' @section Benchmark controls (`bench_config`):
#' * `rounds` (`3L`): number of rounds; odd values let the median ignore a
#'   single cold round.
#' * `iterations` (`1L`): repetitions inside one [bench::mark()] call.
#' * `order` (`"random"`): `"random"` or `"fixed"`.
#' * `seed` (`getFuncOption("seed")`): seed used to shuffle each round.
#' * `warmup` (`TRUE`): run each method once untimed before its first round.
#' * `memory`, `filter_gc`, `check`: forwarded to [bench::mark()].
#' * `auto_budget` (`TRUE`), `mem_safety` (`0.7`), `mem_budget_mb` (`NULL`),
#'   `time_budget_sec` (`NULL`): machine protection.
#' * `measure_rss` (`FALSE`): record peak process RSS (Linux only).
#' * `on_failure` (`"skip"`): `"skip"` or `"stop"`.
#' * `set_threads` (`TRUE`): pin threads to 1 during measurement.
#'
#' @param config A [BenchMemoryConfig] object.
#' @param verbose `logical`. Whether to report progress and print the summary.
#'   Default: `getFuncOption("verbose")`.
#'
#' @return A
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
#'   bench_config = list(rounds = 3L, time_budget_sec = 1800)
#' )
#' res <- BenchMemory(config)
#' attr(res, "summary")
#' }
#'
#' @seealso [BenchMemoryConfig()], [bench::mark()]
#' @family Benchmark
#'
#' @name BenchMemory
#' @export
BenchMemory <- function(obj, ...) {
  UseMethod("BenchMemory")
}

#' @rdname BenchMemory
#' @export
`BenchMemory.SigBridgeR::BenchMemoryConfig` <- function(config, verbose = getFuncOption("verbose")) {
  rlang::check_installed(c("bench", "glue", "callr", "qs2"))

  methods <- unique(config@screen_method)

  if (verbose) {
    cli::cli_progress_step("Detecting machine resources")
  }

  resources <- bench_probe_resources()

  if (verbose) {
    cli::cli_progress_step("Initializing temporary directory")
  }

  # * seed & threads & init
  tmp_dir <- bench_tmp_dir(config@sc_data, config@matched_bulk, config@phenotype)
  on.exit(
    {
      cli::cli_progress_done()
      unlink(tmp_dir, recursive = TRUE, force = TRUE)
      gc(verbose = FALSE)
    },
    add = TRUE
  )

  rows <- vector("list", length(methods) * config@bench_config$rounds)
  rounds <- config@bench_config$rounds

  for (round in seq_len(rounds)) {
    order <- if (config@bench_config$order == "random") {
      set.seed(config@bench_config$seed + round)
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
        config = config,
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
