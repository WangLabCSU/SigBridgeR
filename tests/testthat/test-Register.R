# 80-Register.R ---------------------------------------------------------------

# Registration mutates package-level registry environments in place. Snapshot
# every binding before a test and restore it afterwards so tests stay isolated.
local_registry_snapshot <- function(envs, .local_envir = parent.frame()) {
  before <- lapply(envs, function(env) as.list(env, all.names = TRUE))
  withr::defer(
    {
      for (i in seq_along(envs)) {
        rm(list = ls(envs[[i]], all.names = TRUE), envir = envs[[i]])
        list2env(before[[i]], envir = envs[[i]])
      }
    },
    envir = .local_envir
  )
}

new_screen_method <- function(method_name, method_version = "1.0.0") {
  ScreenMethod(
    method_name = method_name,
    method_version = method_version,
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
}

new_annotation_method <- function(method_name) {
  AnnotationMethod(
    method_name = method_name,
    executor = function(x) NULL,
    pkg_name = "base",
    method_version = "1.0.0"
  )
}

# Register() ------------------------------------------------------------------

test_that("Register() stores ScreenMethod objects under their argument names", {
  local_registry_snapshot(list(ScreenStrategy))
  keys_before <- ls(ScreenStrategy, all.names = TRUE)

  alpha <- new_screen_method("original_alpha")
  beta <- new_screen_method("original_beta")
  Register(alpha = alpha, beta = beta, verbose = FALSE)

  expect_setequal(
    setdiff(ls(ScreenStrategy, all.names = TRUE), keys_before),
    c("alpha", "beta")
  )
  expect_identical(ScreenStrategy$alpha@method_name, "alpha")
  expect_identical(ScreenStrategy$beta@method_name, "beta")
  expect_identical(ScreenStrategy$alpha@executor, alpha@executor)
})

test_that("Register() stores AnnotationMethod objects under their argument names", {
  local_registry_snapshot(list(SCAnnotateStrategy))
  keys_before <- ls(SCAnnotateStrategy, all.names = TRUE)

  ann <- new_annotation_method("original_ann")
  Register(my_ann = ann, verbose = FALSE)

  expect_setequal(
    setdiff(ls(SCAnnotateStrategy, all.names = TRUE), keys_before),
    "my_ann"
  )
  expect_identical(SCAnnotateStrategy$my_ann@method_name, "my_ann")
  expect_identical(SCAnnotateStrategy$my_ann@executor, ann@executor)
})

test_that("Register() stores Seurat functions under single-letter keys", {
  skip_if_not_installed("Seurat")
  local_registry_snapshot(list(SCPreProcessStrategy))
  keys_before <- ls(SCPreProcessStrategy, all.names = TRUE)

  Register(z = Seurat::NormalizeData, verbose = FALSE)

  expect_setequal(
    setdiff(ls(SCPreProcessStrategy, all.names = TRUE), keys_before),
    "z"
  )
  expect_identical(SCPreProcessStrategy$z, Seurat::NormalizeData)
})

test_that("Register() dispatches mixed inputs in a single call", {
  skip_if_not_installed("Seurat")
  local_registry_snapshot(list(
    ScreenStrategy,
    SCAnnotateStrategy,
    SCPreProcessStrategy
  ))

  Register(
    mixed_screen = new_screen_method("screen"),
    mixed_ann = new_annotation_method("ann"),
    y = Seurat::NormalizeData,
    verbose = FALSE
  )

  expect_identical(ScreenStrategy$mixed_screen@method_name, "mixed_screen")
  expect_identical(SCAnnotateStrategy$mixed_ann@method_name, "mixed_ann")
  expect_identical(SCPreProcessStrategy$y, Seurat::NormalizeData)
})

test_that("Register() returns invisibly", {
  local_registry_snapshot(list(ScreenStrategy))

  expect_invisible(Register(invisible_key = new_screen_method("x"), verbose = FALSE))
})

test_that("Register() rejects duplicate keys unless overwrite = TRUE", {
  local_registry_snapshot(list(ScreenStrategy))
  Register(dup = new_screen_method("first"), verbose = FALSE)

  expect_error(
    Register(dup = new_screen_method("second", "2.0.0"), verbose = FALSE),
    "already exists",
    class = "rlang_error"
  )
  expect_identical(ScreenStrategy$dup@method_version, "1.0.0")

  Register(
    dup = new_screen_method("second", "2.0.0"),
    overwrite = TRUE,
    verbose = FALSE
  )
  expect_identical(ScreenStrategy$dup@method_name, "dup")
  expect_identical(ScreenStrategy$dup@method_version, "2.0.0")
})

test_that("Register() allows overwriting annotation and function entries", {
  skip_if_not_installed("Seurat")
  local_registry_snapshot(list(SCAnnotateStrategy, SCPreProcessStrategy))

  Register(dup_ann = new_annotation_method("first"), q = Seurat::NormalizeData, verbose = FALSE)
  Register(
    dup_ann = new_annotation_method("second"),
    q = Seurat::ScaleData,
    overwrite = TRUE,
    verbose = FALSE
  )

  expect_identical(SCAnnotateStrategy$dup_ann@method_name, "dup_ann")
  expect_identical(SCPreProcessStrategy$q, Seurat::ScaleData)
})

test_that("Register() announces registrations only when verbose = TRUE", {
  local_registry_snapshot(list(ScreenStrategy))

  expect_message(
    Register(loud = new_screen_method("loud"), verbose = TRUE),
    "Registered.*ScreenStrategy"
  )
  expect_no_message(
    Register(quiet = new_screen_method("quiet"), verbose = FALSE)
  )
})

# RegisterImpl() ---------------------------------------------------------------

test_that("RegisterImpl() dispatches each supported class to its registry", {
  skip_if_not_installed("Seurat")
  local_registry_snapshot(list(
    ScreenStrategy,
    SCAnnotateStrategy,
    SCPreProcessStrategy
  ))

  screen <- new_screen_method("screen")
  ann <- new_annotation_method("ann")

  expect_message(RegisterImpl(screen, name = "impl_screen", verbose = TRUE), "ScreenStrategy")
  expect_message(RegisterImpl(ann, name = "impl_ann", verbose = TRUE), "SCAnnotateStrategy")
  expect_message(
    RegisterImpl(Seurat::NormalizeData, name = "w", verbose = TRUE),
    "SCPreProcessStrategy"
  )

  expect_identical(ScreenStrategy$impl_screen@method_name, "impl_screen")
  expect_identical(SCAnnotateStrategy$impl_ann@method_name, "impl_ann")
  expect_identical(SCPreProcessStrategy$w, Seurat::NormalizeData)
})

test_that("RegisterImpl() returns invisibly TRUE for supported inputs", {
  local_registry_snapshot(list(ScreenStrategy))

  expect_identical(
    withVisible(RegisterImpl(new_screen_method("ret"), name = "ret", verbose = FALSE)),
    list(value = TRUE, visible = FALSE)
  )
})

test_that("RegisterImpl() rejects unsupported classes", {
  expect_error(
    RegisterImpl(1L, name = "x"),
    "Unsupported class",
    class = "rlang_error"
  )
  expect_error(
    RegisterImpl(list(), name = "x"),
    "[TYPE ERROR]",
    class = "rlang_error",
    fixed = TRUE
  )
})

test_that("RegisterImpl() requires a single-letter key for functions", {
  skip_if_not_installed("Seurat")

  expect_error(
    RegisterImpl(Seurat::NormalizeData, name = "abc", verbose = FALSE),
    "single letter",
    class = "rlang_error"
  )
  expect_error(
    RegisterImpl(Seurat::NormalizeData, name = 1L, verbose = FALSE),
    class = "rlang_error"
  )
})

test_that("RegisterImpl() rejects functions not from Seurat or SeuratObject", {
  expect_error(
    RegisterImpl(function(x) x, name = "z", verbose = FALSE),
    "not a Seurat function",
    class = "rlang_error"
  )
  expect_error(
    RegisterImpl(stats::median, name = "z", verbose = FALSE),
    "not a Seurat function",
    class = "rlang_error"
  )
})

test_that("RegisterImpl() rejects duplicate function keys unless overwrite = TRUE", {
  skip_if_not_installed("Seurat")
  local_registry_snapshot(list(SCPreProcessStrategy))

  RegisterImpl(Seurat::NormalizeData, name = "x", verbose = FALSE)
  expect_error(
    RegisterImpl(Seurat::ScaleData, name = "x", verbose = FALSE),
    "already exists",
    class = "rlang_error"
  )
  expect_identical(SCPreProcessStrategy$x, Seurat::NormalizeData)

  RegisterImpl(Seurat::ScaleData, name = "x", overwrite = TRUE, verbose = FALSE)
  expect_identical(SCPreProcessStrategy$x, Seurat::ScaleData)
})

# is_func_from_pkg() -----------------------------------------------------------

test_that("is_func_from_pkg() detects functions from the requested namespace", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")

  expect_true(is_func_from_pkg(Seurat::NormalizeData, "Seurat"))
  expect_true(is_func_from_pkg(SeuratObject::CreateSeuratObject, "SeuratObject"))
})

test_that("is_func_from_pkg() returns FALSE for local and foreign functions", {
  skip_if_not_installed("Seurat")

  expect_false(is_func_from_pkg(function(x) x, "Seurat"))
  expect_false(is_func_from_pkg(stats::median, "Seurat"))
})

test_that("is_func_from_pkg() validates its inputs", {
  expect_error(is_func_from_pkg(1L, "Seurat"), "must be a function")
  expect_error(
    is_func_from_pkg(function(x) x, 1L),
    "must be character"
  )
})
