test_that("Register works", {
  Register(h=Seurat::Seurat::RunCCA)

  Register(my = DoScissor)

  Register(R=Seurat::RunCCA, my = DoScissor)
})
