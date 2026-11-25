# Tests for the S7 property definitions and abstract base class in
# R/04a-ClassBase.R.
#
# The properties are internal building blocks for SigBridgeRBase subclasses.
# These tests cover:
#   * the declared S7 class and default of every property;
#   * the contract enforced by each property validator;
#   * propagation of S7 class and validator errors through real classes;
#   * the abstractness and version default of SigBridgeRBase.
#
# Implementation notes asserted below:
#   * The list properties surface the last validator condition only.
#     `property_list` therefore rejects duplicate names but accepts empty and
#     unnamed lists; `property_data_list` accepts NULL, NA, empty and unnamed
#     lists. ScreenMethodResult relies on `extra = list()` being accepted.
#   * The version properties are not actually read-only: a valid dotted string
#     can still be assigned after construction.

# ---- helpers ----------------------------------------------------------------

# S7 class names are global, so use a fresh name for every test-file run.
new_test_child <- function(properties = list(), prefix = "TestSigBridgeRBase") {
  S7::new_class(
    name = paste0(prefix, "Child", sample.int(1e9, 1L)),
    parent = SigBridgeRBase,
    properties = properties
  )
}

assign_prop <- function(object, name, value) {
  S7::prop(object, name) <- value
  object
}

fake_seurat <- structure(list(), class = "Seurat")

# A subclass exercising every property defined in this file.
full_child_class <- new_test_child(
  properties = list(
    chr = property_chr,
    phenotype_class = property_phenotype_class,
    fn = property_fn,
    lst = property_list,
    mapper = property_mapper_fn,
    py_env = property_py_env,
    data_list = property_data_list,
    seurat = property_seurat,
    phenotype = property_phenotype,
    matrix = property_matrix,
    version = property_verison
  ),
  prefix = "TestSigBridgeRBaseFull"
)

# ---- version properties -----------------------------------------------------

test_that("property_sigbridger_verison defaults to the package version", {
  expect_identical(property_sigbridger_verison$class, S7::class_character)
  expect_identical(property_sigbridger_verison$default, get_pkg_version())
  expect_type(property_sigbridger_verison$default, "character")
  expect_match(property_sigbridger_verison$default, "\\.")
  expect_null(
    property_sigbridger_verison$validator(
      property_sigbridger_verison$default
    )
  )
})

test_that("property_sigbridger_verison rejects strings without a dot", {
  expect_null(property_sigbridger_verison$validator("4.0"))
  expect_null(property_sigbridger_verison$validator("4.0.0"))

  for (bad in c("4", "", "noversion")) {
    expect_match(
      property_sigbridger_verison$validator(bad),
      "does not contain a dot",
      fixed = TRUE
    )
  }
})

test_that("property_verison is optional but rejects dotless strings", {
  expect_identical(property_verison$class, S7::class_any)
  expect_null(property_verison$default)
  expect_null(property_verison$validator(NULL))
  expect_null(property_verison$validator("4.0.0"))
  expect_null(property_verison$validator("4.0"))

  expect_match(
    property_verison$validator("4"),
    "does not contain a dot",
    fixed = TRUE
  )
  expect_match(
    property_verison$validator(""),
    "does not contain a dot",
    fixed = TRUE
  )
})

# ---- simple typed properties ------------------------------------------------

test_that("plain typed properties have the expected S7 class", {
  expect_identical(property_chr$class, S7::class_character)
  expect_null(property_chr$validator)

  expect_identical(property_fn$class, S7::class_function)
  expect_null(property_fn$validator)
})

test_that("property_phenotype_class accepts only supported labels", {
  expect_identical(property_phenotype_class$class, S7::class_character)
  expect_null(property_phenotype_class$validator("binary"))
  expect_null(property_phenotype_class$validator(c("binary", "survival")))
  expect_null(
    property_phenotype_class$validator(
      c("binary", "continuous", "survival")
    )
  )

  for (bad in list("categorical", c("binary", "unknown"), NA_character_)) {
    expect_match(
      property_phenotype_class$validator(bad),
      "must be one or more of",
      fixed = TRUE
    )
  }
})

test_that("property_mapper_fn accepts functions or NULL", {
  expect_identical(property_mapper_fn$class, S7::class_any)
  expect_null(property_mapper_fn$validator(identity))
  expect_null(property_mapper_fn$validator(\(x) x))
  expect_null(property_mapper_fn$validator(NULL))

  for (bad in list(1L, "identity", list())) {
    expect_match(
      property_mapper_fn$validator(bad),
      "function must be a",
      fixed = TRUE
    )
  }
})

test_that("property_py_env accepts characters or NULL", {
  expect_identical(property_py_env$class, S7::class_any)
  expect_null(property_py_env$validator("my-env"))
  expect_null(property_py_env$validator(c("conda", "path")))
  expect_null(property_py_env$validator(NULL))
  expect_null(property_py_env$validator(NA_character_))

  for (bad in list(42L, TRUE, list("env"))) {
    expect_match(property_py_env$validator(bad), "must be a", fixed = TRUE)
    expect_match(property_py_env$validator(bad), "character", fixed = TRUE)
  }
})

# ---- list properties --------------------------------------------------------

test_that("property_list rejects duplicate names", {
  expect_identical(property_list$class, S7::class_list)
  expect_null(property_list$validator(list(a = 1L, b = 2L)))

  expect_match(
    property_list$validator(list(a = 1L, a = 2L)),
    "names must be unique",
    fixed = TRUE
  )
  expect_match(
    property_list$validator(list(a = 1L, a = 2L, b = 3L)),
    "names must be unique",
    fixed = TRUE
  )
})

test_that("property_list currently accepts empty and unnamed lists", {
  # The naming/missing-value check is not returned by the validator, so only
  # the uniqueness check below it can fail. Keep this behaviour stable because
  # ScreenMethodResult passes `extra = list(...)`, which is `list()` by default.
  expect_null(property_list$validator(list()))
  expect_null(property_list$validator(list(1L, 2L)))
  expect_null(property_list$validator(list(a = 1L, 2L)))
  expect_null(property_list$validator(NULL))
  expect_null(property_list$validator(NA))
})

test_that("property_data_list allows optional and unnamed lists", {
  expect_identical(property_data_list$class, S7::class_list)
  expect_null(property_data_list$validator(NULL))
  expect_null(property_data_list$validator(NA))
  expect_null(property_data_list$validator(list()))
  expect_null(property_data_list$validator(list(1L, 2L)))
  expect_null(property_data_list$validator(list(a = 1L, b = 2L)))

  expect_match(
    property_data_list$validator(list(a = 1L, a = 2L)),
    "names must be unique",
    fixed = TRUE
  )
})

# ---- object-valued properties -----------------------------------------------

test_that("property_seurat requires an object inheriting from Seurat", {
  expect_identical(property_seurat$class, S7::class_any)
  expect_null(property_seurat$validator(fake_seurat))

  skip_if_not_installed("SeuratObject")
  seurat_obj <- suppressWarnings(
    SeuratObject::CreateSeuratObject(counts = matrix(1:4, nrow = 2L))
  )
  expect_null(property_seurat$validator(seurat_obj))

  for (bad in list(NULL, data.frame(x = 1), list())) {
    expect_match(
      property_seurat$validator(bad),
      "must be a",
      fixed = TRUE
    )
    expect_match(property_seurat$validator(bad), "Seurat", fixed = TRUE)
  }
})

test_that("property_phenotype accepts named vectors and 2-D phenotype data", {
  expect_identical(property_phenotype$class, S7::class_any)

  expect_null(property_phenotype$validator(c(time = 1, status = 0)))
  expect_null(property_phenotype$validator(c(a = "x", b = "y")))
  expect_null(
    property_phenotype$validator(data.frame(time = 1, status = 0))
  )
  expect_null(
    property_phenotype$validator(
      data.frame(time = 1, status = 0, score = 2)
    )
  )
  expect_null(
    property_phenotype$validator(
      matrix(
        1:4,
        nrow = 2L,
        dimnames = list(NULL, c("time", "status"))
      )
    )
  )
})

test_that("property_phenotype rejects unnamed vectors and incomplete data", {
  unnamed_vector <- c(1, 2)
  expect_match(
    property_phenotype$validator(unnamed_vector),
    "named vector",
    fixed = TRUE
  )

  incomplete_df <- data.frame(score = 1)
  expect_match(
    property_phenotype$validator(incomplete_df),
    "data.frame",
    fixed = TRUE
  )

  matrix_without_columns <- matrix(1:4, nrow = 2L)
  expect_match(
    property_phenotype$validator(matrix_without_columns),
    "data.frame",
    fixed = TRUE
  )

  expect_match(
    property_phenotype$validator(NULL),
    "named vector",
    fixed = TRUE
  )
})

test_that("property_matrix requires a base matrix", {
  expect_identical(property_matrix$class, S7::class_any)
  expect_null(property_matrix$validator(matrix(1:4, nrow = 2L)))
  expect_null(property_matrix$validator(matrix(letters[1:4], nrow = 2L)))
  expect_null(property_matrix$validator(matrix(TRUE, nrow = 2L, ncol = 2L)))

  for (bad in list(1:4, data.frame(x = 1), NULL)) {
    expect_match(
      property_matrix$validator(bad),
      "must be a",
      fixed = TRUE
    )
    expect_match(property_matrix$validator(bad), "matrix", fixed = TRUE)
  }
})

# ---- SigBridgeRBase ---------------------------------------------------------

test_that("SigBridgeRBase is abstract and declares its version property", {
  expect_s3_class(SigBridgeRBase, "S7_class")
  expect_identical(SigBridgeRBase@name, "SigBridgeRBase")
  expect_true(SigBridgeRBase@abstract)
  expect_identical(
    names(SigBridgeRBase@properties),
    "sigbridger_version"
  )
  version_prop <- SigBridgeRBase@properties$sigbridger_version
  expect_identical(version_prop$class, property_sigbridger_verison$class)
  expect_identical(version_prop$default, property_sigbridger_verison$default)
  expect_identical(version_prop$validator, property_sigbridger_verison$validator)
  expect_error(SigBridgeRBase(), "abstract", fixed = TRUE)
})

test_that("SigBridgeRBase subclasses inherit the package version", {
  child_class <- new_test_child()
  obj <- child_class()

  expect_true(S7::S7_inherits(obj, SigBridgeRBase))
  expect_identical(names(obj), "sigbridger_version")
  expect_identical(obj@sigbridger_version, get_pkg_version())
  expect_identical(
    S7::prop(obj, "sigbridger_version"),
    get_pkg_version()
  )
})

test_that("SigBridgeRBase subclasses validate the version on assignment", {
  child_class <- new_test_child()
  obj <- child_class()

  expect_error(
    assign_prop(obj, "sigbridger_version", "noversion"),
    "does not contain a dot",
    fixed = TRUE
  )

  # The property is mutable despite the documentation calling it read-only.
  obj <- assign_prop(obj, "sigbridger_version", "9.9.9")
  expect_identical(obj@sigbridger_version, "9.9.9")
})

# ---- S7 integration ---------------------------------------------------------

test_that("a full subclass accepts valid values for every property", {
  obj <- full_child_class(
    chr = "value",
    phenotype_class = c("binary", "continuous"),
    fn = identity,
    lst = list(a = 1L),
    mapper = NULL,
    py_env = "sigbridger-env",
    data_list = list(a = 1L),
    seurat = fake_seurat,
    phenotype = c(time = 1, status = 0),
    matrix = matrix(1:4, nrow = 2L),
    version = "1.2.3"
  )

  expect_true(S7::S7_inherits(obj, SigBridgeRBase))
  expect_identical(obj@chr, "value")
  expect_identical(obj@phenotype_class, c("binary", "continuous"))
  expect_identical(obj@fn, identity)
  expect_identical(obj@lst, list(a = 1L))
  expect_null(obj@mapper)
  expect_identical(obj@py_env, "sigbridger-env")
  expect_identical(obj@data_list, list(a = 1L))
  expect_true(inherits(obj@seurat, "Seurat"))
  expect_identical(obj@phenotype, c(time = 1, status = 0))
  expect_identical(obj@matrix, matrix(1:4, nrow = 2L))
  expect_identical(obj@version, "1.2.3")
  expect_identical(obj@sigbridger_version, get_pkg_version())
})

test_that("S7 class checks reject values of the wrong type", {
  obj <- full_child_class(
    chr = "value",
    phenotype_class = "binary",
    fn = identity,
    lst = list(a = 1L),
    mapper = NULL,
    py_env = NULL,
    data_list = list(a = 1L),
    seurat = fake_seurat,
    phenotype = c(time = 1, status = 0),
    matrix = matrix(1:4, nrow = 2L),
    version = "1.2.3"
  )

  expect_error(
    assign_prop(obj, "chr", 1L),
    "must be <character>",
    fixed = TRUE
  )
  expect_error(
    assign_prop(obj, "fn", 1L),
    "must be <function>",
    fixed = TRUE
  )
})

test_that("S7 property validators reject invalid assignments", {
  obj <- full_child_class(
    chr = "value",
    phenotype_class = "binary",
    fn = identity,
    lst = list(a = 1L),
    mapper = NULL,
    py_env = NULL,
    data_list = list(a = 1L),
    seurat = fake_seurat,
    phenotype = c(time = 1, status = 0),
    matrix = matrix(1:4, nrow = 2L),
    version = "1.2.3"
  )

  expect_error(
    assign_prop(obj, "phenotype_class", "categorical"),
    "must be one or more of",
    fixed = TRUE
  )
  expect_error(
    assign_prop(obj, "lst", list(a = 1L, a = 2L)),
    "names must be unique",
    fixed = TRUE
  )
  expect_error(
    assign_prop(obj, "mapper", 1L),
    "function must be a",
    fixed = TRUE
  )
  expect_error(
    assign_prop(obj, "py_env", 1L),
    "must be a <character>",
    fixed = TRUE
  )
  expect_error(
    assign_prop(obj, "data_list", list(a = 1L, a = 2L)),
    "names must be unique",
    fixed = TRUE
  )
  expect_error(
    assign_prop(obj, "seurat", data.frame(x = 1)),
    "must be a <Seurat>",
    fixed = TRUE
  )
  expect_error(
    assign_prop(obj, "phenotype", c(1, 2)),
    "named vector",
    fixed = TRUE
  )
  expect_error(
    assign_prop(obj, "phenotype", data.frame(score = 1)),
    "data.frame",
    fixed = TRUE
  )
  expect_error(
    assign_prop(obj, "matrix", 1:4),
    "must be a <matrix>",
    fixed = TRUE
  )
  expect_error(
    assign_prop(obj, "version", "1"),
    "does not contain a dot",
    fixed = TRUE
  )
})
