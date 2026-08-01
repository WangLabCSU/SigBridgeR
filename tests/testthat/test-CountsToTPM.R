# CountsToTPM - base matrix input ----------------------------------------

describe("CountsToTPM - base matrix input", {
  it("converts simple matrix to TPM (columns sum to 1e6)", {
    counts <- matrix(
      c(100, 200, 300, 150, 250, 350),
      nrow = 3,
      ncol = 2,
      dimnames = list(c("G1", "G2", "G3"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000, G2 = 2000, G3 = 1500)

    result <- CountsToTPM(counts, gene_length)

    expect_equal(dim(result), c(3, 2))
    expect_equal(colSums(result), c(S1 = 1e6, S2 = 1e6), tolerance = 1)
    expect_equal(rownames(result), c("G1", "G2", "G3"))
    expect_equal(colnames(result), c("S1", "S2"))
  })

  it("returns numeric matrix output", {
    counts <- matrix(
      c(100, 200, 150, 250),
      nrow = 2,
      ncol = 2,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000, G2 = 2000)

    result <- CountsToTPM(counts, gene_length)

    expect_true(is.matrix(result))
    expect_type(result, "double")
  })

  it("gene with longer length gets lower TPM for same count", {
    counts <- matrix(
      c(100, 100),
      nrow = 2,
      ncol = 1,
      dimnames = list(c("short", "long"), "S1")
    )
    gene_length <- c(short = 1000, long = 10000)

    result <- CountsToTPM(counts, gene_length)

    # Same raw count, but long gene should have lower TPM (divided by length)
    expect_gt(result["short", "S1"], result["long", "S1"])
  })

  it("handles single gene", {
    counts <- matrix(
      c(100, 200),
      nrow = 1,
      ncol = 2,
      dimnames = list("G1", c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000)

    result <- CountsToTPM(counts, gene_length)

    expect_equal(dim(result), c(1, 2))
    # Each column should sum to 1e6 (only one gene)
    expect_equal(unname(colSums(result)), c(1e6, 1e6), tolerance = 1)
  })

  it("handles single sample", {
    counts <- matrix(
      c(100, 200, 300),
      nrow = 3,
      ncol = 1,
      dimnames = list(c("G1", "G2", "G3"), "S1")
    )
    gene_length <- c(G1 = 1000, G2 = 2000, G3 = 1500)

    result <- CountsToTPM(counts, gene_length)

    expect_equal(dim(result), c(3, 1))
    expect_equal(sum(result), 1e6, tolerance = 1)
  })

  it("handles all-zero column (returns zeros for that column)", {
    counts <- matrix(
      c(0, 0, 100, 200),
      nrow = 2,
      ncol = 2,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000, G2 = 2000)

    result <- CountsToTPM(counts, gene_length)

    # Column S1: all zeros -> zero denominator -> all zeros
    expect_equal(result[, "S1"], c(G1 = 0, G2 = 0))
    # Column S2: should still sum to 1e6
    expect_equal(sum(result[, "S2"]), 1e6, tolerance = 1)
  })

  it("accepts gene_length as a named list", {
    counts <- matrix(
      c(100, 200),
      nrow = 2,
      ncol = 1,
      dimnames = list(c("G1", "G2"), "S1")
    )
    gene_length <- list(G1 = 1000, G2 = 2000)

    result <- CountsToTPM(counts, gene_length)

    expect_equal(sum(result), 1e6, tolerance = 1)
  })

  it("preserves dimnames in output", {
    counts <- matrix(
      c(100, 200, 150, 250),
      nrow = 2,
      ncol = 2,
      dimnames = list(c("GeneA", "GeneB"), c("SampleX", "SampleY"))
    )
    gene_length <- c(GeneA = 1000, GeneB = 2000)

    result <- CountsToTPM(counts, gene_length)

    expect_equal(rownames(result), c("GeneA", "GeneB"))
    expect_equal(colnames(result), c("SampleX", "SampleY"))
  })
})

# CountsToTPM - data.frame input -----------------------------------------

describe("CountsToTPM - data.frame input", {
  it("converts data.frame to matrix and computes TPM", {
    counts <- data.frame(
      S1 = c(100, 200),
      S2 = c(150, 250),
      row.names = c("G1", "G2")
    )
    gene_length <- c(G1 = 1000, G2 = 2000)

    result <- CountsToTPM(counts, gene_length)

    expect_equal(dim(result), c(2, 2))
    expect_equal(unname(colSums(result)), c(1e6, 1e6), tolerance = 1)
  })
})

# CountsToTPM - sparse matrix input (dgCMatrix) --------------------------

describe("CountsToTPM - dgCMatrix input", {
  it("converts dgCMatrix to TPM and returns dgCMatrix", {
    skip_if_not_installed("Matrix")
    counts <- Matrix::Matrix(
      c(100, 0, 200, 0, 300, 0),
      nrow = 3,
      ncol = 2,
      sparse = TRUE,
      dimnames = list(c("G1", "G2", "G3"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000, G2 = 2000, G3 = 1500)

    result <- CountsToTPM(counts, gene_length)

    expect_s4_class(result, "dgCMatrix")
    expect_equal(dim(result), c(3, 2))
  })

  it("dgCMatrix columns sum to 1e6 (non-zero columns)", {
    skip_if_not_installed("Matrix")
    counts <- Matrix::Matrix(
      c(100, 200, 300, 150, 250, 350),
      nrow = 3,
      ncol = 2,
      sparse = TRUE,
      dimnames = list(c("G1", "G2", "G3"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000, G2 = 2000, G3 = 1500)

    result <- CountsToTPM(counts, gene_length)

    expect_equal(unname(Matrix::colSums(result)), c(1e6, 1e6), tolerance = 1)
  })

  it("dgCMatrix preserves sparsity structure", {
    skip_if_not_installed("Matrix")
    counts <- Matrix::Matrix(
      c(100, 0, 200, 0, 300, 0),
      nrow = 3,
      ncol = 2,
      sparse = TRUE,
      dimnames = list(c("G1", "G2", "G3"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000, G2 = 2000, G3 = 1500)

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
      c(100, 200, 150, 250),
      nrow = 2,
      ncol = 2,
      sparse = FALSE,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000, G2 = 2000)

    result <- CountsToTPM(counts, gene_length)

    expect_s4_class(result, "dgeMatrix")
    expect_equal(dim(result), c(2, 2))
    expect_equal(Matrix::colSums(result), c(S1 = 1e6, S2 = 1e6), tolerance = 1)
  })
})

# CountsToTPM - gene_length validation -----------------------------------

describe("CountsToTPM - gene_length validation", {
  it("aborts when gene_length is unnamed", {
    counts <- matrix(
      1:4,
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )

    expect_error(
      CountsToTPM(counts, c(1000, 2000)),
      "must be named"
    )
  })

  it("aborts when gene_length has empty names", {
    counts <- matrix(
      1:4,
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000, 2000) # second has empty name

    expect_error(
      CountsToTPM(counts, gene_length),
      "empty names"
    )
  })

  it("aborts when gene_length has non-positive values", {
    counts <- matrix(
      1:4,
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 0, G2 = 2000)

    expect_error(
      CountsToTPM(counts, gene_length),
      "positive finite"
    )
  })

  it("aborts when gene_length has negative values", {
    counts <- matrix(
      1:4,
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = -100, G2 = 2000)

    expect_error(
      CountsToTPM(counts, gene_length),
      "positive finite"
    )
  })

  it("aborts when gene_length has NA values", {
    counts <- matrix(
      1:4,
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = NA_real_, G2 = 2000)

    expect_error(
      CountsToTPM(counts, gene_length),
      "positive finite"
    )
  })

  it("aborts when gene_length has Inf values", {
    counts <- matrix(
      1:4,
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = Inf, G2 = 2000)

    expect_error(
      CountsToTPM(counts, gene_length),
      "positive finite"
    )
  })

  it("aborts when gene_length missing a gene in counts", {
    counts <- matrix(
      1:4,
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000) # G2 missing

    expect_error(
      CountsToTPM(counts, gene_length),
      "misses gene"
    )
  })
})

# CountsToTPM - counts validation ----------------------------------------

describe("CountsToTPM - counts validation", {
  it("aborts when counts is not matrix/dgCMatrix/dgeMatrix", {
    gene_length <- c(G1 = 1000, G2 = 2000)

    expect_error(
      CountsToTPM(c(100, 200), gene_length),
      "Expected.*counts"
    )
  })

  it("aborts when counts has no rownames", {
    counts <- matrix(1:4, nrow = 2, ncol = 2)
    gene_length <- c(G1 = 1000, G2 = 2000)

    expect_error(
      CountsToTPM(counts, gene_length),
      "rownames"
    )
  })

  it("aborts when counts contains NA values", {
    counts <- matrix(
      c(100, NA, 150, 250),
      nrow = 2,
      ncol = 2,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000, G2 = 2000)

    expect_error(
      CountsToTPM(counts, gene_length),
      "NA"
    )
  })

  it("aborts when counts contains negative values", {
    counts <- matrix(
      c(100, -5, 150, 250),
      nrow = 2,
      ncol = 2,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000, G2 = 2000)

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
      c(50, 150, 250, 80, 200, 320),
      nrow = 3,
      ncol = 2,
      dimnames = list(c("G1", "G2", "G3"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 500, G2 = 1500, G3 = 2500)

    result <- CountsToTPM(counts, gene_length)

    expect_equal(colSums(result), c(S1 = 1e6, S2 = 1e6), tolerance = 1)
  })

  it("TPM values are proportional to counts/length", {
    counts <- matrix(
      c(100, 100, 200, 200),
      nrow = 2,
      ncol = 2,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
    gene_length <- c(G1 = 1000, G2 = 1000)

    result <- CountsToTPM(counts, gene_length)

    # Same length, count ratio 1:2 -> TPM ratio 1:2
    expect_equal(result["G2", "S1"] / result["G1", "S1"], 1)
    expect_equal(result["G2", "S2"] / result["G1", "S2"], 1)
  })

  it("larger counts matrix still works correctly", {
    set.seed(42)
    n_genes <- 100
    n_samples <- 5
    counts <- matrix(
      rpois(n_genes * n_samples, lambda = 100),
      nrow = n_genes,
      ncol = n_samples,
      dimnames = list(
        paste0("G", seq_len(n_genes)),
        paste0("S", seq_len(n_samples))
      )
    )
    gene_length <- setNames(
      runif(n_genes, 500, 5000),
      paste0("G", seq_len(n_genes))
    )

    result <- CountsToTPM(counts, gene_length)

    expect_equal(dim(result), c(n_genes, n_samples))
    expect_equal(unname(colSums(result)), rep(1e6, n_samples), tolerance = 1)
  })
})
