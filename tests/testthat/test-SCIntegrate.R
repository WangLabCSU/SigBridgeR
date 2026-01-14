test_that("It works with matrices", {
  mat1 <- matrix(
    rpois(100, 5),
    nrow = 20,
    dimnames = list(paste0("G", 1:20), paste0("C", 1:5))
  )
  mat2 <- matrix(
    rpois(120, 6),
    nrow = 20,
    dimnames = list(paste0("G", 11:30), paste0("C", 1:6))
  )

  # * check name generating
  quo <- rlang::enquos(mat1, mat2)
  expect_vector(get_names_4_ids(quo))

  integrated <- SCIntegrate(mat1, mat2)
  # colnames(integrated)

  expect_type(integrated, "integer")
  expect_equal(ncol(integrated), 11)
  expect_equal(nrow(integrated), 30)

  integrated <- SCIntegrate(A = mat1, B = mat2)
  # colnames(integrated)

  mat3 <- Matrix::Matrix(mat1)
  mat4 <- Matrix::Matrix(mat2)

  integrated2 <- SCIntegrate(mat3, mat4)
  class(integrated2)
  expect_s4_class(integrated2, "dgCMatrix")
})

test_that("It works with Seurat", {
  mat1 <- matrix(
    rpois(1000, 5),
    nrow = 20,
    dimnames = list(paste0("G", 1:20), paste0("C", 1:50))
  )
  mat2 <- matrix(
    rpois(1200, 6),
    nrow = 20,
    dimnames = list(paste0("G", 11:30), paste0("C", 1:60))
  )

  seu1 <- Seurat::CreateSeuratObject(mat1)
  seu2 <- Seurat::CreateSeuratObject(mat2)
  integrated_seu <- SCIntegrate(
    seu1,
    seu2,
    method = Seurat::CCAIntegration,
    pipeline = "nsvpi",
    dims = 1:10,
    k.weight = 40
  )

  expect_s4_class(integrated_seu, "Seurat")
})
