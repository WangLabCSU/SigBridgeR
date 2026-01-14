test_that("multiplication works", {
  bulk <- readRDS("/data/resource/wanglab/TCGA/count/TCGA-CESC.exp.count.rds")

  head(bulk)

  filtered <- SigBridgeR::BulkPreProcess(bulk)

  head(filtered)

  transformed <- log2(filtered + 1)

  head(transformed)
})
