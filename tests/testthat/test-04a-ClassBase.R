test_that("version properties validate dotted version strings", {
  expect_identical(property_sigbridger_verison$class, S7::class_character)
  expect_identical(property_sigbridger_verison$default, get_pkg_version())
  expect_null(property_sigbridger_verison$validator(get_pkg_version()))

  expect_identical(property_verison$class, S7::class_any)
  expect_null(property_verison$default)
  expect_null(property_verison$validator(NULL))
  expect_null(property_verison$validator("4.0.0"))
  expect_match(
    property_verison$validator("4"),
    "does not contain a dot",
    fixed = TRUE
  )
})

test_that("plain typed properties have the expected S7 class", {
  expect_identical(property_chr$class, S7::class_character)
  expect_null(property_chr$validator)

  expect_identical(property_fn$class, S7::class_function)
  expect_null(property_fn$validator)
})

test_that("property_phenotype_class accepts only supported phenotype classes", {
  expect_identical(property_phenotype_class$class, S7::class_character)
  expect_null(property_phenotype_class$validator("binary"))
  expect_null(property_phenotype_class$validator(c("binary", "survival")))

  expect_match(
    property_phenotype_class$validator("categorical"),
    "must be one or more of",
    fixed = TRUE
  )
  expect_match(
    property_phenotype_class$validator(c("binary", "unknown")),
    "must be one or more of",
    fixed = TRUE
  )
})

test_that("property_list validates name uniqueness", {
  expect_identical(property_list$class, S7::class_list)
  expect_null(property_list$validator(list(a = 1L, b = 2L)))
  expect_null(property_list$validator(list()))
  expect_null(property_list$validator(list(1L, 2L)))
  expect_null(property_list$validator(list(a = 1L, 2L)))

  expect_match(
    property_list$validator(list(a = 1L, a = 2L)),
    "names must be unique",
    fixed = TRUE
  )
})

test_that("property_mapper_fn accepts functions or NULL", {
  expect_identical(property_mapper_fn$class, S7::class_any)
  expect_null(property_mapper_fn$validator(function(x) x))
  expect_null(property_mapper_fn$validator(NULL))

  msg <- property_mapper_fn$validator(1L)
  expect_match(msg, "function must be a", fixed = TRUE)
  expect_match(msg, "function", fixed = TRUE)
})

test_that("property_py_env accepts characters or NULL", {
  expect_identical(property_py_env$class, S7::class_any)
  expect_null(property_py_env$validator("my-env"))
  expect_null(property_py_env$validator(NULL))

  msg <- property_py_env$validator(42L)
  expect_match(msg, "must be a", fixed = TRUE)
  expect_match(msg, "character", fixed = TRUE)
})

test_that("property_data_list requires unique names but allows loose names", {
  expect_identical(property_data_list$class, S7::class_list)
  expect_null(property_data_list$validator(list(a = 1L, b = 2L)))
  expect_null(property_data_list$validator(list()))
  expect_null(property_data_list$validator(list(1L, 2L)))

  expect_match(
    property_data_list$validator(list(a = 1L, a = 2L)),
    "names must be unique",
    fixed = TRUE
  )
})

test_that("property_seurat requires a Seurat object", {
  expect_identical(property_seurat$class, S7::class_any)
  expect_null(property_seurat$validator(structure(list(), class = "Seurat")))

  msg <- property_seurat$validator(data.frame(x = 1))
  expect_match(msg, "value must be a", fixed = TRUE)
  expect_match(msg, "Seurat", fixed = TRUE)
})

test_that("property_phenotype accepts named vectors and two-dimensional time/status data", {
  expect_identical(property_phenotype$class, S7::class_any)
  expect_null(property_phenotype$validator(c(time = 1, status = 0)))
  expect_null(property_phenotype$validator(data.frame(time = 1, status = 0)))
  expect_null(
    property_phenotype$validator(
      matrix(
        1:4,
        nrow = 2L,
        dimnames = list(NULL, c("time", "status"))
      )
    )
  )

  msg_vector <- property_phenotype$validator(c(1, 2))
  expect_match(msg_vector, "named vector", fixed = TRUE)

  msg_df <- property_phenotype$validator(data.frame(score = 1))
  expect_match(msg_df, "data.frame", fixed = TRUE)

  msg_matrix <- property_phenotype$validator(matrix(1:4, nrow = 2L))
  expect_match(msg_matrix, "data.frame", fixed = TRUE)
})

test_that("property_matrix requires a matrix", {
  expect_identical(property_matrix$class, S7::class_any)
  expect_null(property_matrix$validator(matrix(1:4, nrow = 2L)))

  msg <- property_matrix$validator(1:4)
  expect_match(msg, "value must be a", fixed = TRUE)
  expect_match(msg, "matrix", fixed = TRUE)
})

test_that("SigBridgeRBase is abstract and provides a validated version property", {
  expect_identical(SigBridgeRBase@name, "SigBridgeRBase")
  expect_true(SigBridgeRBase@abstract)
  expect_identical(names(SigBridgeRBase@properties), "sigbridger_version")
  expect_error(SigBridgeRBase(), "abstract")

  child_name <- paste0("TestSigBridgeRBaseChild", sample.int(1e6, 1L))
  child_class <- S7::new_class(child_name, parent = SigBridgeRBase)
  obj <- child_class()

  expect_true(S7::S7_inherits(obj, SigBridgeRBase))
  expect_identical(obj@sigbridger_version, get_pkg_version())
  expect_error(
    S7::prop(obj, "sigbridger_version") <- "noversion",
    "does not contain a dot",
    fixed = TRUE
  )
})
