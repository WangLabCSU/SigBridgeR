test_that("DrawBulkPCA returns invisible plot object", {
  set.seed(123)

  n_genes <- 50L
  n_samples <- 30L
  bulk <- matrix(
    rnorm(n_genes * n_samples, mean = 10, sd = 2),
    nrow = n_genes,
    ncol = n_samples,
    dimnames = list(
      paste0("Gene", seq_len(n_genes)),
      paste0("Sample", seq_len(n_samples))
    )
  )

  group <- rep(c("Control", "Treatment"), times = c(15L, 15L))

  # show_plot = FALSE should return invisibly
  p <- DrawBulkPCA(bulk = bulk, group = group, show_plot = FALSE)
  expect_s3_class(p, "patchwork")
})

test_that("DrawBulkPCA works with batch argument", {
  set.seed(123)

  n_genes <- 20L
  n_samples <- 20L
  bulk <- matrix(
    rnorm(n_genes * n_samples, mean = 10, sd = 2),
    nrow = n_genes,
    ncol = n_samples,
    dimnames = list(
      paste0("Gene", seq_len(n_genes)),
      paste0("Sample", seq_len(n_samples))
    )
  )

  group <- rep(c("A", "B"), each = 10L)
  batch <- rep(c("Batch1", "Batch2", "Batch3"), length.out = n_samples)

  p <- DrawBulkPCA(
    bulk = bulk,
    group = group,
    batch = batch,
    show_plot = FALSE
  )
  expect_s3_class(p, "patchwork")
})

test_that("DrawBulkPCA errors when batch has too many levels", {
  set.seed(123)

  n_genes <- 10L
  n_samples <- 10L
  bulk <- matrix(
    rnorm(n_genes * n_samples),
    nrow = n_genes,
    ncol = n_samples,
    dimnames = list(
      paste0("Gene", seq_len(n_genes)),
      paste0("Sample", seq_len(n_samples))
    )
  )

  group <- rep("A", n_samples)
  batch <- paste0("Batch", seq_len(n_samples))

  expect_error(
    DrawBulkPCA(bulk = bulk, group = group, batch = batch, show_plot = FALSE),
    class = "chk_error"
  )
})
