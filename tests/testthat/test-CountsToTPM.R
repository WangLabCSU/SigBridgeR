# CountsToTPM - base matrix input ----------------------------------------

describe("CountsToTPM - base matrix input", {
  it("converts simple matrix to TPM (columns sum to 1e6)", {
    counts <- matrix(
      c(100L, 200L, 300L, 150L, 250L, 350L),
      nrow = 3L,
      ncol = 2L,
      dimnames = list(c("G1", "G2", "G3"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000L, G2 = 2000L, G3 = 1500L)

    result <- CountsToTPM(counts, gene_length)

    expect_equal(dim(result), c(3L, 2L))
    expect_equal(colSums(result), c(S1 = 1e6, S2 = 1e6), tolerance = 1L)
    expect_equal(rownames(result), c("G1", "G2", "G3"))
    expect_equal(colnames(result), c("S1", "S2"))
  })

  it("returns numeric matrix output", {
    counts <- matrix(
      c(100L, 200L, 150L, 250L),
      nrow = 2L,
      ncol = 2L,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000L, G2 = 2000L)

    result <- CountsToTPM(counts, gene_length)

    expect_true(is.matrix(result))
    expect_type(result, "double")
  })

  it("gene with longer length gets lower TPM for same count", {
    counts <- matrix(
      c(100L, 100L),
      nrow = 2L,
      ncol = 1L,
      dimnames = list(c("short", "long"), "S1")
    )
    gene_length <- c(short = 1000L, long = 10000L)

    result <- CountsToTPM(counts, gene_length)

    # Same raw count, but long gene should have lower TPM (divided by length)
    expect_gt(result["short", "S1"], result["long", "S1"])
  })

  it("handles single gene", {
    counts <- matrix(
      c(100L, 200L),
      nrow = 1L,
      ncol = 2L,
      dimnames = list("G1", c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000L)

    result <- CountsToTPM(counts, gene_length)

    expect_equal(dim(result), c(1L, 2L))
    # Each column should sum to 1e6 (only one gene)
    expect_equal(unname(colSums(result)), c(1e6, 1e6), tolerance = 1L)
  })

  it("handles single sample", {
    counts <- matrix(
      c(100L, 200L, 300L),
      nrow = 3L,
      ncol = 1L,
      dimnames = list(c("G1", "G2", "G3"), "S1")
    )
    gene_length <- c(G1 = 1000L, G2 = 2000L, G3 = 1500L)

    result <- CountsToTPM(counts, gene_length)

    expect_equal(dim(result), c(3L, 1L))
    expect_equal(sum(result), 1e6, tolerance = 1L)
  })

  it("handles all-zero column (returns zeros for that column)", {
    counts <- matrix(
      c(0L, 0L, 100L, 200L),
      nrow = 2L,
      ncol = 2L,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000L, G2 = 2000L)

    result <- CountsToTPM(counts, gene_length)

    # Column S1: all zeros -> zero denominator -> all zeros
    expect_equal(result[, "S1"], c(G1 = 0L, G2 = 0L))
    # Column S2: should still sum to 1e6
    expect_equal(sum(result[, "S2"]), 1e6, tolerance = 1L)
  })

  it("accepts gene_length as a named list", {
    counts <- matrix(
      c(100L, 200L),
      nrow = 2L,
      ncol = 1L,
      dimnames = list(c("G1", "G2"), "S1")
    )
    gene_length <- list(G1 = 1000L, G2 = 2000L)

    result <- CountsToTPM(counts, gene_length)

    expect_equal(sum(result), 1e6, tolerance = 1L)
  })

  it("preserves dimnames in output", {
    counts <- matrix(
      c(100L, 200L, 150L, 250L),
      nrow = 2L,
      ncol = 2L,
      dimnames = list(c("GeneA", "GeneB"), c("SampleX", "SampleY"))
    )
    gene_length <- c(GeneA = 1000L, GeneB = 2000L)

    result <- CountsToTPM(counts, gene_length)

    expect_equal(rownames(result), c("GeneA", "GeneB"))
    expect_equal(colnames(result), c("SampleX", "SampleY"))
  })
})

# CountsToTPM - data.frame input -----------------------------------------

describe("CountsToTPM - data.frame input", {
  it("converts data.frame to matrix and computes TPM", {
    counts <- data.frame(
      S1 = c(100L, 200L),
      S2 = c(150L, 250L),
      row.names = c("G1", "G2")
    )
    gene_length <- c(G1 = 1000L, G2 = 2000L)

    result <- CountsToTPM(counts, gene_length)

    expect_equal(dim(result), c(2L, 2L))
    expect_equal(unname(colSums(result)), c(1e6, 1e6), tolerance = 1L)
  })
})

# CountsToTPM - sparse matrix input (dgCMatrix) --------------------------

describe("CountsToTPM - dgCMatrix input", {
  it("converts dgCMatrix to TPM and returns dgCMatrix", {
    skip_if_not_installed("Matrix")
    counts <- Matrix::Matrix(
      c(100L, 0L, 200L, 0L, 300L, 0L),
      nrow = 3L,
      ncol = 2L,
      sparse = TRUE,
      dimnames = list(c("G1", "G2", "G3"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000L, G2 = 2000L, G3 = 1500L)

    result <- CountsToTPM(counts, gene_length)

    expect_s4_class(result, "Matrix")
    expect_equal(dim(result), c(3L, 2L))
  })

  it("dgCMatrix columns sum to 1e6 (non-zero columns)", {
    skip_if_not_installed("Matrix")
    counts <- Matrix::Matrix(
      c(100L, 200L, 300L, 150L, 250L, 350L),
      nrow = 3L,
      ncol = 2L,
      sparse = TRUE,
      dimnames = list(c("G1", "G2", "G3"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000L, G2 = 2000L, G3 = 1500L)

    result <- CountsToTPM(counts, gene_length)

    expect_equal(unname(Matrix::colSums(result)), c(1e6, 1e6), tolerance = 1L)
  })

  it("dgCMatrix preserves sparsity structure", {
    skip_if_not_installed("Matrix")
    counts <- Matrix::Matrix(
      c(100L, 0L, 200L, 0L, 300L, 0L),
      nrow = 3L,
      ncol = 2L,
      sparse = TRUE,
      dimnames = list(c("G1", "G2", "G3"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000L, G2 = 2000L, G3 = 1500L)

    result <- CountsToTPM(counts, gene_length)

    # Same sparsity pattern
    expect_equal(Matrix::nnzero(result), Matrix::nnzero(counts))
  })
})

# CountsToTPM - dgeMatrix input ------------------------------------------

describe("CountsToTPM - dgeMatrix input", {
  it("converts dgeMatrix to TPM", {
    skip_if_not_installed("Matrix")
    counts <- Matrix::Matrix(
      c(100L, 200L, 150L, 250L),
      nrow = 2L,
      ncol = 2L,
      sparse = FALSE,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000L, G2 = 2000L)

    result <- CountsToTPM(counts, gene_length)

    expect_s4_class(result, "dgeMatrix")
    expect_equal(dim(result), c(2L, 2L))
    expect_equal(Matrix::colSums(result), c(S1 = 1e6, S2 = 1e6), tolerance = 1L)
  })
})

# CountsToTPM - gene_length validation -----------------------------------

describe("CountsToTPM - gene_length validation", {
  it("aborts when gene_length is unnamed", {
    counts <- matrix(
      1L:4L,
      nrow = 2L,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )

    expect_error(
      CountsToTPM(counts, c(1000L, 2000L)),
      "must be a named vector"
    )
  })

  it("aborts when gene_length has empty names", {
    counts <- matrix(
      1L:4L,
      nrow = 2L,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000L, 2000L) # second has empty name

    expect_error(
      CountsToTPM(counts, gene_length),
      "empty names"
    )
  })

  it("aborts when gene_length has non-positive values", {
    counts <- matrix(
      1L:4L,
      nrow = 2L,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 0L, G2 = 2000L)

    expect_error(
      CountsToTPM(counts, gene_length),
      "positive finite"
    )
  })

  it("aborts when gene_length has negative values", {
    counts <- matrix(
      1L:4L,
      nrow = 2L,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = -100L, G2 = 2000L)

    expect_error(
      CountsToTPM(counts, gene_length),
      "positive finite"
    )
  })

  it("aborts when gene_length has NA values", {
    counts <- matrix(
      1L:4L,
      nrow = 2L,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = NA_real_, G2 = 2000L)

    expect_error(
      CountsToTPM(counts, gene_length),
      "positive finite"
    )
  })

  it("aborts when gene_length has Inf values", {
    counts <- matrix(
      1L:4L,
      nrow = 2L,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = Inf, G2 = 2000L)

    expect_error(
      CountsToTPM(counts, gene_length),
      "positive finite"
    )
  })

  it("aborts when gene_length missing a gene in counts", {
    counts <- matrix(
      1L:4L,
      nrow = 2L,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000L) # G2 missing

    expect_error(
      CountsToTPM(counts, gene_length),
      "misses gene"
    )
  })
})

# CountsToTPM - counts validation ----------------------------------------

describe("CountsToTPM - counts validation", {
  it("aborts when counts is not matrix/dgCMatrix/dgeMatrix", {
    gene_length <- c(G1 = 1000L, G2 = 2000L)

    expect_error(
      CountsToTPM(c(100L, 200L), gene_length),
      "must have rownames"
    )
  })

  it("aborts when counts has no rownames", {
    counts <- matrix(1L:4L, nrow = 2L, ncol = 2L)
    gene_length <- c(G1 = 1000L, G2 = 2000L)

    expect_error(
      CountsToTPM(counts, gene_length),
      "rownames"
    )
  })

  it("aborts when counts contains NA values", {
    counts <- matrix(
      c(100L, NA, 150L, 250L),
      nrow = 2L,
      ncol = 2L,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000L, G2 = 2000L)

    expect_error(
      CountsToTPM(counts, gene_length),
      "NA"
    )
  })

  it("aborts when counts contains negative values", {
    counts <- matrix(
      c(100L, -5L, 150L, 250L),
      nrow = 2L,
      ncol = 2L,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000L, G2 = 2000L)

    expect_error(
      CountsToTPM(counts, gene_length),
      "non-negative"
    )
  })
})

# CountsToTPM - TPM properties -------------------------------------------

describe("CountsToTPM - TPM properties", {
  it("each column sums to 1e6", {
    counts <- matrix(
      c(50L, 150L, 250L, 80L, 200L, 320L),
      nrow = 3L,
      ncol = 2L,
      dimnames = list(c("G1", "G2", "G3"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 500L, G2 = 1500L, G3 = 2500L)

    result <- CountsToTPM(counts, gene_length)

    expect_equal(colSums(result), c(S1 = 1e6, S2 = 1e6), tolerance = 1L)
  })

  it("TPM values are proportional to counts/length", {
    counts <- matrix(
      c(100L, 100L, 200L, 200L),
      nrow = 2L,
      ncol = 2L,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000L, G2 = 1000L)

    result <- CountsToTPM(counts, gene_length)

    # Same length, count ratio 1:2 -> TPM ratio 1:2
    expect_equal(result["G2", "S1"] / result["G1", "S1"], 1L)
    expect_equal(result["G2", "S2"] / result["G1", "S2"], 1L)
  })

  it("larger counts matrix still works correctly", {
    set.seed(42L)
    n_genes <- 100L
    n_samples <- 5L
    counts <- matrix(
      rpois(n_genes * n_samples, lambda = 100L),
      nrow = n_genes,
      ncol = n_samples,
      dimnames = list(
        paste0("G", seq_len(n_genes)),
        paste0("S", seq_len(n_samples))
      )
    )
    gene_length <- setNames(
      runif(n_genes, 500L, 5000L),
      paste0("G", seq_len(n_genes))
    )

    result <- CountsToTPM(counts, gene_length)

    expect_equal(dim(result), c(n_genes, n_samples))
    expect_equal(unname(colSums(result)), rep(1e6, n_samples), tolerance = 1L)
  })
})
