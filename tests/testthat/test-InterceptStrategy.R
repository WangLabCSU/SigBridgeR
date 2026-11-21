# 80-InterceptStrategy.R ----------------------------------------------------

test_that("InterceptStrategy() summarizes an environment of plain functions", {
  skip_if_not_installed("dplyr")

  env <- rlang::new_environment(list(
    alpha = function(x, y = 1) NULL,
    beta = function(z) NULL
  ))

  info <- InterceptStrategy(env)

  expect_s3_class(info, "tbl_df")
  expect_named(info, c("name", "args"))
  expect_setequal(info$name, c("alpha", "beta"))

  args <- setNames(info$args, info$name)
  expect_identical(args$alpha, c("x", "y"))
  expect_identical(args$beta, "z")
})

test_that("InterceptStrategy() resolves full and partially matched strategy names", {
  skip_if_not_installed("dplyr")

  full <- InterceptStrategy("SCPreProcessStrategy")
  partial <- InterceptStrategy("SCPreProcess")

  expect_identical(full, partial)
  expect_named(full, c("name", "args"))
  expect_setequal(full$name, names(SCPreProcessStrategy))
  expect_type(full$args, "list")
  expect_true(all(vapply(full$args, is.character, logical(1))))
})

test_that("InterceptStrategy() inspects all three registered strategy environments", {
  skip_if_not_installed("dplyr")

  screen <- InterceptStrategy("ScreenStrategy")
  preprocess <- InterceptStrategy("SCPreProcessStrategy")
  annotate <- InterceptStrategy("SCAnnotateStrategy")

  expect_setequal(screen$name, AVAILABLE_SCREEN_METHODS)
  expect_setequal(preprocess$name, names(SCPreProcessStrategy))
  expect_setequal(annotate$name, AVAILABLE_ANNOTATE_METHODS)

  expect_true(all(c("executor", "phenotype_class") %in% names(screen)))
  expect_false("phenotype_class" %in% names(annotate))
  expect_true(all(c("name", "args") %in% names(preprocess)))
})

test_that("InterceptStrategy() expands S7 method properties into list columns", {
  skip_if_not_installed("dplyr")

  method <- ScreenMethod(
    method_name = "TestMethod",
    method_version = "1.0.0",
    executor = function(
      matched_bulk,
      sc_data,
      phenotype,
      label_type,
      phenotype_class
    ) {
      NULL
    },
    phenotype_class = "binary",
    mapper = NULL
  )
  info <- InterceptStrategy(rlang::new_environment(list(test = method)))

  expect_s3_class(info, "tbl_df")
  expect_identical(info$name, "test")
  expect_identical(info$method_name, "TestMethod")
  expect_identical(
    info$executor[[1L]],
    c("matched_bulk", "sc_data", "phenotype", "label_type", "phenotype_class")
  )
  expect_identical(info$phenotype_class[[1L]], "binary")
})

test_that("InterceptStrategy() aborts when the strategy has no names", {
  skip_if_not_installed("dplyr")

  expect_error(
    InterceptStrategy(list(function(x) NULL)),
    class = "rlang_error"
  )
  expect_error(
    InterceptStrategy(list(function(x) NULL)),
    "CONFIG ERROR",
    class = "rlang_error"
  )
})

test_that("InterceptStrategy() rejects additional arguments", {
  skip_if_not_installed("dplyr")

  expect_error(
    InterceptStrategy(
      rlang::new_environment(list(f = function(x) NULL)),
      extra = TRUE
    ),
    "must be empty",
    class = "rlang_error"
  )
})

test_that("get_strategy_info() records the name and formals of a function", {
  info <- get_strategy_info(function(a, b = 2) NULL, "fn")

  expect_named(info, c("name", "args"))
  expect_identical(info$name, "fn")
  expect_identical(info$args, list(c("a", "b")))
})

test_that("get_strategy_info() wraps executor formals and phenotype_class for S7 methods", {
  method <- ScreenMethod(
    method_name = "TestMethod",
    method_version = "1.0.0",
    executor = function(
      matched_bulk,
      sc_data,
      phenotype,
      label_type,
      phenotype_class
    ) {
      NULL
    },
    phenotype_class = c("binary", "survival"),
    mapper = NULL
  )

  info <- get_strategy_info(method, "method")

  expect_identical(
    info$executor,
    list(
      c(
        "matched_bulk",
        "sc_data",
        "phenotype",
        "label_type",
        "phenotype_class"
      )
    )
  )
  expect_identical(info$phenotype_class, list(c("binary", "survival")))
  expect_identical(info$name, "method")
  expect_identical(info$method_name, "TestMethod")
})

test_that("get_strategy_info() skips missing executor or phenotype_class properties", {
  config <- ScreenMethodConfig(
    method_name = "Scissor",
    param = list(alpha = 0.05),
    phenotype_class = "binary",
    label_type = "binary",
    method_version = "1.0"
  )
  config_info <- get_strategy_info(config, "config")

  expect_null(config_info$executor)
  expect_identical(config_info$phenotype_class, list("binary"))
  expect_identical(config_info$name, "config")

  annotation <- AnnotationMethod(
    method_name = "A",
    executor = function(x) NULL,
    method_version = "1.0"
  )
  annotation_info <- get_strategy_info(annotation, "annotation")

  expect_identical(annotation_info$executor, list("x"))
  expect_null(annotation_info$phenotype_class)
  expect_identical(annotation_info$name, "annotation")
})
