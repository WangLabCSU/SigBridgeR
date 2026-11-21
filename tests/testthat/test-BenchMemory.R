test_that("BenchMemoryConfigSet merges defaults and rejects unknown fields", {
  ctrl <- BenchMemoryConfigSet()
  expect_equal(ctrl$rounds, 3L)
  expect_equal(ctrl$order, "random")
  expect_true(ctrl$memory)
  expect_equal(ctrl$timeout_sec, 10800L)

  ctrl <- BenchMemoryConfigSet(rounds = 5L, order = "fixed")
  expect_equal(ctrl$rounds, 5L)
  expect_equal(ctrl$order, "fixed")

  expect_error(BenchMemoryConfigSet(other_xyz = 5L), "Unknown")
  expect_error(BenchMemoryConfigSet(rounds = 0L), "at least 1")
})

# Mocked measurement: the orchestration is tested without launching child R sessions.
mock_measure <- function(seen) {
  function(config, method_name, round, position, tmp_dir, verbose = FALSE) {
    seen$value <- c(seen$value, paste0(round, ":", method_name))
    data.table::data.table(
      method_name = method_name,
      round = as.integer(round),
      position = as.integer(position),
      status = if (identical(method_name, "broken")) "error" else "ok",
      error = if (identical(method_name, "broken")) "boom" else NA_character_,
      mem_peak_mb = 1,
      rss_peak_mb = NA_real_,
      elapsed = 0.1
    )
  }
}

test_that("BenchMemoryConfig keeps label templates and per-method arguments", {
  skip_if_not_installed("SeuratObject")

  sc <- suppressWarnings(
    SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(matrix(1:4, 2L)))
  )
  config <- BenchMemoryConfig(
    matched_bulk = matrix(1:4, nrow = 2L, dimnames = list(c("g1", "g2"), c("s1", "s2"))),
    sc_data = sc,
    phenotype = c(s1 = 0, s2 = 1),
    phenotype_class = "binary",
    label_type = "{method_name}_res",
    screen_method = "scPP",
    method_param = MethodParamSet(scPP_param = list(Log2FC_cutoff = 0.05)),
    bench_config = BenchMemoryConfigSet(rounds = 1L)
  )

  expect_identical(config@label_type, "{method_name}_res")
  expect_equal(config@method_param[["scPP_param"]]$Log2FC_cutoff, 0.05)
  expect_identical(config@bench_config$order, "random")
})

test_that("BenchMemory runs each method once per round in fixed order", {
  skip_if_not_installed("bench")
  skip_if_not_installed("SeuratObject")

  fake_executor <- function(
    matched_bulk,
    sc_data,
    phenotype,
    label_type,
    phenotype_class,
    ...
  ) {
    list(scRNA_data = sc_data)
  }
  for (nm in c("cheap", "broken")) {
    ScreenStrategy[[nm]] <- ScreenMethod(
      method_name = nm,
      method_version = "0.0.1",
      executor = fake_executor,
      phenotype_class = "binary",
      mapper = NULL
    )
  }
  withr::defer({
    rm("cheap", envir = ScreenStrategy)
    rm("broken", envir = ScreenStrategy)
  })

  sc <- suppressWarnings(
    SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(matrix(1:4, 2L)))
  )
  config <- BenchMemoryConfig(
    matched_bulk = matrix(1:4, nrow = 2L, dimnames = list(c("g1", "g2"), c("s1", "s2"))),
    sc_data = sc,
    phenotype = c(s1 = 0, s2 = 1),
    phenotype_class = "binary",
    screen_method = c("cheap", "broken"),
    method_param = MethodParamSet(),
    bench_config = BenchMemoryConfigSet(rounds = 2L, order = "fixed")
  )

  seen <- new.env(parent = emptyenv())
  seen$value <- character()
  local_mocked_bindings(
    bench_probe_resources = function() list(cores = 1L),
    bench_tmp_dir = function(sc_data, bulk_data, pheno_data) tempfile(),
    bench_mem_measure = mock_measure(seen),
    .package = "SigBridgeR"
  )

  res <- BenchMemory(config, verbose = FALSE)

  expect_s3_class(res, "data.table")
  expect_equal(nrow(res), 4L)
  expect_identical(seen$value, c("1:cheap", "1:broken", "2:cheap", "2:broken"))
  expect_setequal(res$status, c("ok", "error"))
  expect_true(all(res[method_name == "broken" & status == "error"]$error == "boom"))
})

test_that("BenchMemory random order is reproducible for a fixed seed", {
  skip_if_not_installed("bench")
  skip_if_not_installed("SeuratObject")

  fake_executor <- function(
    matched_bulk,
    sc_data,
    phenotype,
    label_type,
    phenotype_class,
    ...
  ) {
    list(scRNA_data = sc_data)
  }
  methods <- c("alpha", "beta", "gamma")
  for (nm in methods) {
    ScreenStrategy[[nm]] <- ScreenMethod(
      method_name = nm,
      method_version = "0.0.1",
      executor = fake_executor,
      phenotype_class = "binary",
      mapper = NULL
    )
  }
  withr::defer(
    for (nm in methods) {
      rm(list = nm, envir = ScreenStrategy)
    }
  )

  sc <- suppressWarnings(
    SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(matrix(1:4, 2L)))
  )

  run_order <- function(seed) {
    config <- BenchMemoryConfig(
      matched_bulk = matrix(1:4, nrow = 2L, dimnames = list(c("g1", "g2"), c("s1", "s2"))),
      sc_data = sc,
      phenotype = c(s1 = 0, s2 = 1),
      phenotype_class = "binary",
      screen_method = methods,
      method_param = MethodParamSet(),
      bench_config = BenchMemoryConfigSet(rounds = 2L, order = "random", seed = seed)
    )
    seen <- new.env(parent = emptyenv())
    seen$value <- character()
    local_mocked_bindings(
      bench_probe_resources = function() list(cores = 1L),
      bench_tmp_dir = function(sc_data, bulk_data, pheno_data) tempfile(),
      bench_mem_measure = mock_measure(seen),
      .package = "SigBridgeR"
    )
    invisible(BenchMemory(config, verbose = FALSE))
    seen$value
  }

  first <- run_order(42L)
  second <- run_order(42L)
  expect_identical(first, second)
  expect_setequal(first[grepl("^1:", first)], paste0("1:", methods))
  expect_setequal(first[grepl("^2:", first)], paste0("2:", methods))
})

test_that("BenchMemory.default reports unsupported classes", {
  expect_error(BenchMemory(list()), "not implemented for class")
})
