test_that("BenchCtrl merges defaults and rejects unknown fields", {
  ctrl <- BenchCtrl()
  expect_equal(ctrl$rounds, 3L)
  expect_equal(ctrl$order, "random")
  expect_true(ctrl$warmup)

  ctrl <- BenchCtrl(list(rounds = 5, order = "fixed"))
  expect_equal(ctrl$rounds, 5L)
  expect_equal(ctrl$order, "fixed")

  expect_error(BenchCtrl(list(round = 5)), "Unknown")
  expect_error(BenchCtrl(list(rounds = 0L)), "at least 1")
})

test_that("bench_round_order is reproducible and honours fixed order", {
  methods <- c("Scissor", "scPP", "scPAS", "LP_SGL")

  expect_identical(bench_round_order(methods, "fixed", 1L), methods)

  first <- bench_round_order(methods, "random", 42L)
  second <- bench_round_order(methods, "random", 42L)
  expect_identical(first, second)
  expect_setequal(first, methods)

  expect_identical(bench_round_order("scPP", "random", 1L), "scPP")
})

test_that("bench_resolve_label interpolates the method name", {
  expect_equal(bench_resolve_label(NULL, "scPP"), "scPP")
  expect_equal(bench_resolve_label("{method_name}", "scPP"), "scPP")
  expect_equal(bench_resolve_label("{method_name}_res", "scPAS"), "scPAS_res")
  expect_equal(bench_resolve_label("label", "scPAS"), "label")
})

test_that("bench_summarise aggregates rounds per method", {
  metrics <- data.table::data.table(
    method_name = c("A", "A", "A", "B"),
    method_version = "1.0",
    round = c(1L, 2L, 3L, 1L),
    status = "ok",
    time_sec = c(1, 2, 3, 10),
    mem_alloc_mb = c(1, 1, 1, 2),
    mem_peak_mb = c(10, 12, 11, 20),
    rss_peak_mb = NA_real_
  )

  res <- bench_summarise(metrics)
  expect_equal(nrow(res), 2L)
  expect_equal(res[method_name == "A"]$time_median_sec, 2)
  expect_equal(res[method_name == "A"]$time_max_sec, 3)
  expect_equal(res[method_name == "A"]$n_rounds, 3L)
  expect_true(res[1L]$method_name == "A")
})

test_that("bench_bind_metrics keeps failed runs visible", {
  rows <- list(data.table::data.table(
    method_name = "A",
    round = 1L,
    status = "ok",
    time_sec = 1
  ))
  errors <- list(data.table::data.table(
    method_name = "B",
    round = 1L,
    status = "error",
    error = "boom"
  ))

  res <- bench_bind_metrics(rows, errors)
  expect_equal(nrow(res), 2L)
  expect_setequal(res$method_name, c("A", "B"))
  expect_equal(res[method_name == "B"]$error, "boom")
})

test_that("bench_method_params splices method specific arguments", {
  fake_executor <- function(
    matched_bulk,
    sc_data,
    phenotype,
    label_type,
    phenotype_class,
    ...
  ) {
    list(scRNA_data = sc_data, ...)
  }

  ScreenStrategy[["fake"]] <- ScreenMethod(
    method_name = "fake",
    method_version = "0.0.1",
    executor = fake_executor,
    phenotype_class = c("binary", "survival"),
    mapper = function(params) {
      params$mapped <- TRUE
      params
    }
  )
  withr::defer(rm("fake", envir = ScreenStrategy))

  skip_if_not_installed("SeuratObject")
  sc <- SeuratObject::CreateSeuratObject(
    counts = matrix(1:4, nrow = 2, dimnames = list(c("g1", "g2"), c("c1", "c2")))
  )

  config <- BenchMemoryConfig(
    matched_bulk = matrix(1:4, nrow = 2, dimnames = list(c("g1", "g2"), c("s1", "s2"))),
    sc_data = sc,
    phenotype = c(s1 = 0, s2 = 1),
    screen_method = "fake",
    method_param = MethodParamSet(list(fake_param = list(alpha = 0.05)))
  )

  params <- bench_method_params(config, "fake")
  expect_equal(params$label_type, "fake")
  expect_equal(params$phenotype_class, "binary")
  expect_equal(params$alpha, 0.05)
  expect_true(params$mapped)
})

test_that("BenchMemory measures every method and returns a bench object", {
  skip_if_not_installed("bench")
  skip_if_not_installed("SeuratObject")

  cheap_executor <- function(
    matched_bulk,
    sc_data,
    phenotype,
    label_type,
    phenotype_class,
    ...
  ) {
    list(scRNA_data = sc_data)
  }
  failing_executor <- function(
    matched_bulk,
    sc_data,
    phenotype,
    label_type,
    phenotype_class,
    ...
  ) {
    stop("intentional failure")
  }

  ScreenStrategy[["cheap"]] <- ScreenMethod(
    method_name = "cheap",
    method_version = "0.0.1",
    executor = cheap_executor,
    phenotype_class = c("binary", "survival", "continuous"),
    mapper = NULL
  )
  ScreenStrategy[["broken"]] <- ScreenMethod(
    method_name = "broken",
    method_version = "0.0.1",
    executor = failing_executor,
    phenotype_class = c("binary", "survival", "continuous"),
    mapper = NULL
  )
  withr::defer({
    rm("cheap", envir = ScreenStrategy)
    rm("broken", envir = ScreenStrategy)
  })

  sc <- SeuratObject::CreateSeuratObject(
    counts = matrix(1:4, nrow = 2, dimnames = list(c("g1", "g2"), c("c1", "c2")))
  )

  config <- BenchMemoryConfig(
    matched_bulk = matrix(1:4, nrow = 2, dimnames = list(c("g1", "g2"), c("s1", "s2"))),
    sc_data = sc,
    phenotype = c(s1 = 0, s2 = 1),
    screen_method = c("cheap", "broken"),
    method_param = MethodParamSet(),
    bench_config = list(rounds = 2L, warmup = FALSE, auto_budget = FALSE)
  )

  res <- BenchMemory(config, verbose = FALSE)

  expect_s3_class(res, "bench_mark")
  expect_true(all(c("method_name", "round", "time_sec", "mem_peak_mb") %in% names(res)))
  expect_equal(nrow(res[method_name == "cheap"]), 2L)
  expect_equal(nrow(res[method_name == "broken" & status == "error"]), 2L)

  summary <- attr(res, "summary")
  expect_true("cheap" %chin% summary$method_name)
  expect_false("broken" %chin% summary$method_name)

  settings <- attr(res, "settings")
  expect_equal(settings$control$rounds, 2L)
  expect_true(data.table::is.data.table(settings$history))
})
