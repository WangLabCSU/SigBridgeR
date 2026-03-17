test_that("multiplication works", {
   seurat <- qs::qread("/home/data/sigbridger/merged_without_pipet.qs")

   upset <- ScreenUpset(seurat,n_intersection = 63)
})
