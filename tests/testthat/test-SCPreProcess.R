test_that("SCPreProcess works", {
  # 1. 设置随机种子以保证结果可重复
  set.seed(42)

  # 2. 定义数据规模
  num_genes <- 1000
  num_cells <- 200

  data_matrix <- matrix(
    rpois(num_genes * num_cells, lambda = 100),
    nrow = num_genes,
    ncol = num_cells
  )

  meta_info <- data.frame(
    group = rep(c("Control", "Treatment"), each = num_cells / 2),
    batch = sample(c("Batch_A", "Batch_B"), num_cells, replace = TRUE),
    row.names = colnames(data_matrix)
  )

  rownames(data_matrix) <- paste0("GENE", 1:num_genes)
  colnames(data_matrix) <- paste0("CELL", 1:num_cells)

  data_matrix[1:4, 1:4]
  #   CELL1 CELL2 CELL3 CELL4
  # GENE1   113   106    84    99
  # GENE2    94    84    86   110
  # GENE3   100    98   100    98
  # GENE4    88    97    87   106

  seu1 <- Seurat::CreateSeuratObject(counts = data_matrix)

  #   params <- list(o = list(meta.data = meta_info))

  seu2 <- rlang::exec(
    SeuratObject::CreateSeuratObject,
    counts = data_matrix,
    !!!params$o
  )

  #   counts <- Matrix::Matrix(data_matrix)

  seu <- SCPreProcess(sc = data_matrix)
})
