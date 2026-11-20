# Tests for AddMetaFeature() (R/11d-AddMetaFeature.R) ------------------------

new_test_seurat <- function() {
  counts <- matrix(
    c(1L, 2L, 3L, 4L, 5L, 6L),
    nrow = 3L,
    dimnames = list(c("g1", "g2", "g3"), c("s1", "s2"))
  )
  suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts))
}

feature_named <- function(values = c(1, 2, 3)) {
  stats::setNames(values, c("g1", "g2", "g3"))
}

test_that("AddMetaFeature rejects a non-Seurat input", {
  expect_error(AddMetaFeature(list()), class = "chk_error")
  expect_error(AddMetaFeature(data.frame(x = 1)), class = "chk_error")
})

test_that("AddMetaFeature warns and returns the object when no metadata is given", {
  obj <- new_test_seurat()

  expect_warning(res <- AddMetaFeature(obj), "No inputs were provided")
  expect_identical(res, obj)
  expect_visible(suppressWarnings(AddMetaFeature(obj)))
})

test_that("AddMetaFeature requires a dot name for a vector", {
  obj <- new_test_seurat()
  expect_error(AddMetaFeature(obj, c(1, 2, 3)), "Column name is empty")
})

test_that("AddMetaFeature names an unnamed vector with the feature names", {
  obj <- new_test_seurat()

  expect_warning(res <- AddMetaFeature(obj, score = c(10, 20, 30)), "not named")
  expect_equal(res@assays$RNA@meta.data$score, c(10, 20, 30))
})

test_that("AddMetaFeature accepts a feature-named vector without warning", {
  obj <- new_test_seurat()

  expect_no_warning(res <- AddMetaFeature(obj, score = feature_named(c(4, 5, 6))))
  expect_equal(res@assays$RNA@meta.data$score, c(4, 5, 6))
})

test_that("AddMetaFeature rejects a mis-sized unnamed vector", {
  obj <- new_test_seurat()
  expect_error(
    suppressWarnings(AddMetaFeature(obj, score = c(1, 2))),
    "length does not match"
  )
})

test_that("AddMetaFeature adds every column of a data.frame", {
  obj <- new_test_seurat()
  df <- data.frame(
    score1 = c(1, 2, 3),
    score2 = c("a", "b", "c"),
    row.names = c("g1", "g2", "g3")
  )

  expect_no_warning(res <- AddMetaFeature(obj, df))
  md <- res@assays$RNA@meta.data
  expect_setequal(names(md), c("score1", "score2"))
  expect_equal(md$score1, c(1, 2, 3))
  expect_equal(md$score2, c("a", "b", "c"))
})

test_that("AddMetaFeature processes multiple metadata inputs in one call", {
  obj <- new_test_seurat()
  res <- AddMetaFeature(
    obj,
    score1 = feature_named(c(1, 2, 3)),
    score2 = feature_named(c(4, 5, 6))
  )

  expect_setequal(names(res@assays$RNA@meta.data), c("score1", "score2"))
  expect_equal(res@assays$RNA@meta.data$score2, c(4, 5, 6))
})

test_that("AddMetaFeature falls back to the default assay for an unknown assay", {
  obj <- new_test_seurat()
  res <- AddMetaFeature(obj, score = feature_named(), assay = "NOT_AN_ASSAY")

  expect_equal(res@assays$RNA@meta.data$score, c(1, 2, 3))
  expect_false("NOT_AN_ASSAY" %in% names(res@assays))
})

test_that("AddMetaFeature writes to the explicitly requested assay", {
  obj <- new_test_seurat()
  res <- AddMetaFeature(obj, score = feature_named(), assay = "RNA")

  expect_equal(res@assays$RNA@meta.data$score, c(1, 2, 3))
})

test_that("AddMetaFeature returns the modified object invisibly", {
  obj <- new_test_seurat()
  expect_invisible(AddMetaFeature(obj, score = feature_named()))
  expect_length(obj@assays$RNA@meta.data, 0L)
})
