test_that("generate a mock bulk count matrix", {
  set.seed(123)

  # 参数设置
  n_genes <- 10000 # 基因数量
  n_samples_per_batch <- 15 # 每个批次的样本数
  n_batches <- 2
  n_samples <- n_samples_per_batch * n_batches

  # 生成 Ensembl ID
  ens_ids <- paste0("ENSG", sprintf("%011d", 1:n_genes))

  # 生成样本名（带批次信息）
  sample_names <- c(
    paste0("batch1_sample", 1:n_samples_per_batch),
    paste0("batch2_sample", 1:n_samples_per_batch)
  )

  # 1. 生成基础表达强度（log-normal 分布模拟真实计数分布）
  base_means <- exp(rnorm(n_genes, mean = log(10), sd = 1)) # 平均表达水平

  # 2. 为每个样本生成计数：使用负二项分布模拟计数数据（更真实）
  #    负二项分布参数：size（离散度），mu（均值）
  size_param <- 10 # 控制离散程度，越小越分散、零越多

  # 初始化计数矩阵
  counts <- matrix(0L, nrow = n_genes, ncol = n_samples)

  # 3. 引入批次效应：
  #    - batch2 整体 scale 上调 1.5 倍
  #    - 随机 10% 基因在 batch2 中额外上调 3 倍（模拟批次特异高表达）
  batch_effect_multiplier <- rep(1.0, n_samples)
  batch_effect_multiplier[(n_samples_per_batch + 1):n_samples] <- 1.5

  # 随机选择部分基因受更强批次效应影响
  strong_batch_genes <- sample(n_genes, size = round(0.1 * n_genes))
  strong_effect <- matrix(1.0, nrow = n_genes, ncol = n_samples)
  strong_effect[strong_batch_genes, (n_samples_per_batch + 1):n_samples] <- 3.0

  # 4. 逐样本生成计数
  for (j in seq_len(n_samples)) {
    mu <- base_means * batch_effect_multiplier[j] * strong_effect[, j]
    # 使用 rnbinom 生成负二项分布计数
    counts[, j] <- rnbinom(n = n_genes, size = size_param, mu = mu)
  }

  # 5. 强化零膨胀（可选）：随机将部分低表达值设为 0
  zero_inflate_prob <- pmin(0.8 * exp(-counts / 5), 0.9) # 表达越低，越可能被置零
  zero_mask <- matrix(
    runif(n_genes * n_samples) < zero_inflate_prob,
    nrow = n_genes,
    ncol = n_samples
  )
  counts[zero_mask] <- 0L

  # 转换为整数矩阵（确保是 integer 类型，符合 DESeq2 输入要求）
  counts <- as.matrix(counts)
  storage.mode(counts) <- "integer"

  # 设置行名和列名
  rownames(counts) <- ens_ids
  colnames(counts) <- sample_names

  data <- BulkPreProcess(
    counts,
    sample_info = data.frame(
      batch = sample_(c("batch1", "batch2"), n_samples, replace = TRUE),
      sample = sample_names,
      condition = c(
        rep("control", n_samples_per_batch),
        rep("treatment", n_samples_per_batch)
      )
    ),
    gene_symbol_conversion = TRUE,
    check = TRUE,
    min_count_threshold = 0L,
    min_gene_expressed = 0L,
    min_total_reads = 1e4L,
    min_genes_detected = 100L,
    min_correlation = 0.8,
    n_top_genes = 500L,
    show_plot_results = TRUE,
    verbose = TRUE
  )
})
