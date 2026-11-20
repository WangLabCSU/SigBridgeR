# Shared fixtures for screen-method tests.
# All expression matrices are negative-binomial draws with ~100 rows x ~100 cols
# so the suites stay fast.

nb_matrix <- function(nr = 100L, nc = 100L, seed = 123L, prefix_rows = "G", prefix_cols = "S") {
  withr::with_seed(seed, {
    m <- matrix(
      as.integer(stats::rnbinom(nr * nc, size = 5, mu = 20)),
      nrow = nr,
      ncol = nc
    )
  })
  dimnames(m) <- list(
    paste0(prefix_rows, seq_len(nr)),
    paste0(prefix_cols, seq_len(nc))
  )
  m
}

nb_bulk <- function(n_genes = 100L, n_samples = 100L, seed = 123L) {
  nb_matrix(n_genes, n_samples, seed, "G", "S")
}

nb_sc_matrix <- function(n_genes = 100L, n_cells = 100L, seed = 456L) {
  nb_matrix(n_genes, n_cells, seed, "G", "C")
}

nb_seurat <- function(n_genes = 100L, n_cells = 100L, seed = 456L) {
  counts <- Matrix::Matrix(nb_sc_matrix(n_genes, n_cells, seed), sparse = TRUE)
  suppressWarnings(SeuratObject::CreateSeuratObject(
    counts = counts,
    min.cells = 0L,
    min.features = 0L
  ))
}

nb_binary_phenotype <- function(n = 100L, seed = 789L) {
  withr::with_seed(seed, {
    setNames(sample(c(0L, 1L), n, replace = TRUE), paste0("S", seq_len(n)))
  })
}

nb_continuous_phenotype <- function(n = 100L, seed = 789L) {
  withr::with_seed(seed, {
    setNames(round(stats::rnorm(n, mean = 50, sd = 10), 2L), paste0("S", seq_len(n)))
  })
}

nb_survival_phenotype <- function(n = 100L, seed = 789L) {
  withr::with_seed(seed, {
    data.frame(
      time = round(stats::rexp(n, rate = 0.05) + 1, 2L),
      status = sample(c(0L, 1L), n, replace = TRUE),
      row.names = paste0("S", seq_len(n))
    )
  })
}
