test_that("multiplication works", {
  seurat <- qs::qread("/home/data/sigbridger/scab_result.qs")
  seurat <- seurat$scRNA_data

  counts <- SeuratObject::LayerData(seurat, "counts")
  counts[1:4, 1:4]

  seurat1 <- SeuratObject::CreateSeuratObject(counts)

  res <- ChooseNormalization(
    log_norm = Seurat::NormalizeData(seurat1) %>%
      Seurat::ScaleData() %>%
      Seurat::FindVariableFeatures(),
    sctransform = Seurat::SCTransform(seurat1)
  )
})
