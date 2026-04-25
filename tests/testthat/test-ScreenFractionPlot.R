test_that("multiplication works", {
  seurat <- qs::qread("/home/data/sigbridger/merged_without_pipet.qs")

  fraction <- ScreenFractionPlot(seurat, group_by = "seurat_clusters")
})
