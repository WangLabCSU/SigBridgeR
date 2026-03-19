test_that("DrawPCA works", {
  # 设置随机种子以保证可重复性
  set.seed(123)

  # 1. 生成样本数据
  n_samples <- 60
  n_conditions <- 3
  conditions <- rep(c("Control", "Treatment_A", "Treatment_B"), each = 20)

  # 2. 模拟PCA坐标数据（不同条件有不同分布）
  PC1 <- c(
    rnorm(20, mean = -2, sd = 1), # Control组
    rnorm(20, mean = 2, sd = 1), # Treatment_A组
    rnorm(20, mean = 0, sd = 1.5) # Treatment_B组
  )

  PC2 <- c(
    rnorm(20, mean = 1, sd = 0.8),
    rnorm(20, mean = -1, sd = 0.8),
    rnorm(20, mean = 0, sd = 1.2)
  )

  # 3. 创建批次变量（可选）
  batch <- rep(c("Batch_1", "Batch_2", "Batch_3"), times = 20)

  # 4. 构建pca_df数据框
  pca_df <- data.frame(
    PC1 = PC1,
    PC2 = PC2,
    condition = conditions,
    batch = batch,
    stringsAsFactors = FALSE
  )

  # 5. 创建方差解释比例标签
  var_labels <- c(45.32, 28.67) # PC1和PC2的方差解释百分比

  #   7.
  #   测试函数
  DrawPCA(
    pca_data = pca_df,
    var_labels = var_labels,
    show_plot = TRUE
  )

  # 查看数据结构
  str(pca_df)
  summary(pca_df)

  # 查看条件分布
  table(pca_df$condition)
  table(pca_df$batch)
})
