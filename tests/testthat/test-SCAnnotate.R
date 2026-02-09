test_that("SCAnnotate works", {
  seurat <- qs::qread("/home/data/sigbridger/scab_result.qs", nthreads = 4L)
  seurat <- seurat$scRNA_data

  seurat2 <- SCAnnotate(sc = seurat, method = "SingleR")

  seurat2[[]] %>% colnames()

  table_(seurat2$SingleR_labels)
  overview(seurat2$SingleR_delta_next)
  table_(seurat2$SingleR_pruned_labels)
})

test_that("SCAnnotate works", {
  seurat <- qs::qread("/home/data/sigbridger/scab_result.qs", nthreads = 4L)
  seurat <- seurat$scRNA_data

  seurat3 <- SCAnnotate(
    sc = seurat,
    method = "CellTypist",
    model = "Immune_All_Low.pkl",
    download = FALSE,
    celltypist_tools = "/home/yyx/R/Project/R_code/SigBridgeR/inst/python/73-CellTypistAnnotate.py"
  )
})

test_that("SCAnnotate works", {
  seurat <- qs::qread("/home/data/sigbridger/scab_result.qs", nthreads = 4L)
  seurat <- seurat$scRNA_data

  seurat4 <- SCAnnotate(
    sc = seurat,
    method = "mLLMCelltype",
    model = c("deepseek-chat"),
    api_keys = list(
      deepseek = "sk-02cd71d0b3bd40479e5a6c41e84854d5"
    )
  )
  colnames(seurat4[[]])
})
