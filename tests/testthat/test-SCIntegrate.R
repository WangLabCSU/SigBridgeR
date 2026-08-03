# Helpers -----------------------------------------------------------------

# Create a simple numeric matrix with gene-like rownames and cell-like colnames
new_test_matrix <- function(n_genes = 10L, n_cells = 5L, seed = 42L) {
  withr::local_seed(seed)
  matrix(
    rpois(n_genes * n_cells, 5),
    nrow = n_genes,
    dimnames = list(
      paste0("Gene", seq_len(n_genes)),
      paste0("Cell", seq_len(n_cells))
    )
  )
}

# Create a minimal Seurat object
new_test_seurat <- function(n_cells = 10L, n_genes = 10L, seed = 42L) {
  m <- new_test_matrix(n_genes, n_cells, seed)
  SeuratObject::CreateSeuratObject(counts = m)
}

# SCIntegrate - dispatch ---------------------------------------------------

describe("SCIntegrate - dispatch", {
  it("aborts when no arguments are provided", {
    expect_error(
      SCIntegrate(),
      "No arguments provided"
    )
  })

  it("aborts for unsupported class input", {
    expect_error(
      SCIntegrate(list(a = 1)),
      "No implementation for class"
    )
  })
})

# SCIntegrate.matrix -------------------------------------------------------

describe("SCIntegrate.matrix", {
  it("integrates two base matrices with union of genes", {
    mat1 <- new_test_matrix(10L, 5L)
    mat2 <- new_test_matrix(10L, 6L)
    # Shift gene names so there's partial overlap
    rownames(mat2) <- paste0("Gene", 6:15)

    integrated <- SCIntegrate(mat1, mat2)

    expect_s4_class(integrated, "dgCMatrix")
    expect_equal(ncol(integrated), 11L)
    # Union of genes: Gene1-5 (mat1 only) + Gene6-10 (both) + Gene11-15 (mat2 only)
    expect_equal(nrow(integrated), 15L)
  })

  it("returns dgCMatrix for base matrix inputs", {
    mat1 <- new_test_matrix(5L, 3L)
    mat2 <- new_test_matrix(5L, 2L)

    integrated <- SCIntegrate(mat1, mat2)

    expect_s4_class(integrated, "dgCMatrix")
  })

  it("returns dgCMatrix for Matrix (dgCMatrix) inputs", {
    mat1 <- Matrix::Matrix(new_test_matrix(5L, 3L))
    mat2 <- Matrix::Matrix(new_test_matrix(5L, 2L))

    integrated <- SCIntegrate(mat1, mat2)

    expect_s4_class(integrated, "dgCMatrix")
  })

  it("handles mixed base matrix and Matrix inputs", {
    mat1 <- new_test_matrix(5L, 3L)
    mat2 <- Matrix::Matrix(new_test_matrix(5L, 2L))

    integrated <- SCIntegrate(mat1, mat2)

    expect_s4_class(integrated, "dgCMatrix")
    expect_equal(ncol(integrated), 5L)
  })

  it("handles dense Matrix (dgeMatrix) inputs via dispatch", {
    mat1 <- Matrix::Matrix(new_test_matrix(5L, 3L))
    mat2 <- Matrix::Matrix(new_test_matrix(5L, 2L))

    integrated <- SCIntegrate(mat1, mat2)

    expect_s4_class(integrated, "dgCMatrix")
    expect_equal(ncol(integrated), 5L)
    expect_equal(nrow(integrated), 5L)
  })

  it("handles single matrix input", {
    mat <- new_test_matrix(5L, 3L)

    integrated <- SCIntegrate(mat)

    expect_s4_class(integrated, "dgCMatrix")
    expect_equal(ncol(integrated), 3L)
    expect_equal(nrow(integrated), 5L)
  })

  it("respects named arguments as cell prefixes", {
    mat1 <- new_test_matrix(5L, 3L)
    mat2 <- new_test_matrix(5L, 2L)

    integrated <- SCIntegrate(SampleA = mat1, SampleB = mat2)

    expect_true(all(grepl("^SampleA_", colnames(integrated)[1:3])))
    expect_true(all(grepl("^SampleB_", colnames(integrated)[4:5])))
  })

  it("infers prefixes from variable names when unnamed", {
    mat1 <- new_test_matrix(5L, 3L)
    mat2 <- new_test_matrix(5L, 2L)

    integrated <- SCIntegrate(mat1, mat2)

    expect_true(all(grepl("^mat1_", colnames(integrated)[1:3])))
    expect_true(all(grepl("^mat2_", colnames(integrated)[4:5])))
  })

  it("fills missing gene values with 0 (sparse Matrix)", {
    mat1 <- new_test_matrix(3L, 2L)
    rownames(mat1) <- c("A", "B", "C")
    mat2 <- new_test_matrix(3L, 2L)
    rownames(mat2) <- c("B", "C", "D")

    integrated <- SCIntegrate(mat1, mat2)

    # Gene A only in mat1 → present in cols 1-2, 0 in cols 3-4
    expect_equal(as.numeric(integrated["A", 1:2]), as.numeric(mat1["A", ]))
    expect_equal(as.numeric(integrated["A", 3:4]), c(0, 0))
    # Gene D only in mat2 → 0 in cols 1-2, present in cols 3-4
    expect_equal(as.numeric(integrated["D", 1:2]), c(0, 0))
    expect_equal(as.numeric(integrated["D", 3:4]), as.numeric(mat2["D", ]))
  })

  it("preserves data integrity after integration", {
    mat1 <- matrix(
      c(1, 2, 3, 4),
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("C1", "C2"))
    )
    mat2 <- matrix(
      c(5, 6, 7, 8),
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("C3", "C4"))
    )

    integrated <- SCIntegrate(mat1, mat2)

    # colnames get prefixed with argument name
    expect_equal(as.numeric(integrated["G1", "mat1_C1"]), 1)
    expect_equal(as.numeric(integrated["G2", "mat1_C2"]), 4)
    expect_equal(as.numeric(integrated["G1", "mat2_C3"]), 5)
    expect_equal(as.numeric(integrated["G2", "mat2_C4"]), 8)
  })
})

# SCIntegrate.data.frame --------------------------------------------------

describe("SCIntegrate.data.frame", {
  it("integrates two data.frame inputs", {
    # data.frame: rows = genes (rownames), columns = samples
    df1 <- data.frame(
      SampleA = 1:3,
      SampleB = 4:6,
      row.names = c("Gene1", "Gene2", "Gene3")
    )
    df2 <- data.frame(
      SampleC = 7:9,
      SampleD = 10:12,
      row.names = c("Gene1", "Gene3", "Gene4")
    )

    integrated <- SCIntegrate(df1, df2)

    expect_true(is.matrix(integrated))
    # Union of genes: Gene1, Gene2, Gene3, Gene4
    expect_equal(nrow(integrated), 4L)
    expect_equal(ncol(integrated), 4L)
    expect_setequal(rownames(integrated), c("Gene1", "Gene2", "Gene3", "Gene4"))
  })

  it("returns a matrix with NA-filled missing genes", {
    # gene A only in df1, gene B only in df2
    df1 <- data.frame(S1 = 1:2, row.names = c("A", "C"))
    df2 <- data.frame(S2 = 3:4, row.names = c("B", "C"))

    integrated <- SCIntegrate(df1, df2)

    # Gene C (shared) has values in both
    expect_false(anyNA(integrated["C", ]))
    # Gene A (df1 only) has NA in df2 columns
    expect_true(anyNA(integrated["A", ]))
  })

  it("respects named arguments as prefixes in colnames", {
    df1 <- data.frame(S1 = 1:2, row.names = c("X", "Y"))
    df2 <- data.frame(S2 = 3:4, row.names = c("X", "Y"))

    integrated <- SCIntegrate(First = df1, Second = df2)

    expect_true(any(grepl("^First_", colnames(integrated))))
    expect_true(any(grepl("^Second_", colnames(integrated))))
  })
})

# SCIntegrate.Seurat - validation ------------------------------------------

describe("SCIntegrate.Seurat - validation", {
  it("aborts with unknown pipeline steps", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")

    seu1 <- new_test_seurat(10L, 10L)
    seu2 <- new_test_seurat(10L, 10L)

    expect_error(
      SCIntegrate(seu1, seu2, pipeline = "nsvxyz"),
      "Unknown pipeline steps"
    )
  })

  it("warns when method is not a function and returns merged Seurat", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")

    seu1 <- new_test_seurat(10L, 10L)
    seu2 <- new_test_seurat(10L, 10L)

    expect_warning(
      result <- SCIntegrate(seu1, seu2, method = "not_a_function"),
      "method.*must be a function"
    )
    expect_s4_class(result, "Seurat")
  })

  it("aborts when pipeline is not a single character", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")

    seu1 <- new_test_seurat(5L, 5L)
    seu2 <- new_test_seurat(5L, 5L)

    expect_error(
      SCIntegrate(seu1, seu2, pipeline = 123),
    )
  })
})

# SCIntegrate.Seurat - integration -----------------------------------------

describe("SCIntegrate.Seurat - integration", {
  it("integrates two Seurat objects with minimal pipeline and CCAIntegration", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")

    seu1 <- new_test_seurat(50L, 10L)
    seu2 <- new_test_seurat(50L, 10L, seed = 123L)

    integrated <- SCIntegrate(
      seu1,
      seu2,
      method = Seurat::CCAIntegration,
      pipeline = "nsvpi",
      dims = 1:5,
      k.weight = 30
    )

    expect_s4_class(integrated, "Seurat")
  })

  it("respects custom add.cell.ids", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")

    seu1 <- new_test_seurat(5L, 5L)
    seu2 <- new_test_seurat(5L, 5L, seed = 99L)

    # Use a non-function method to skip integration pipeline
    suppressWarnings(
      merged <- SCIntegrate(
        seu1, seu2,
        method = "skip",
        add.cell.ids = c("Batch1", "Batch2")
      )
    )

    expect_true(all(grepl("^Batch1_", colnames(merged)[1:5])))
    expect_true(all(grepl("^Batch2_", colnames(merged)[6:10])))
  })

  it("handles merge.data and merge.dr parameters", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")

    seu1 <- new_test_seurat(5L, 5L)
    seu2 <- new_test_seurat(5L, 5L, seed = 77L)

    suppressWarnings(
      merged <- SCIntegrate(
        seu1, seu2,
        method = "skip",
        merge.data = FALSE,
        merge.dr = NA,
        project = "TestProject"
      )
    )

    expect_s4_class(merged, "Seurat")
    expect_equal(merged@project.name, "TestProject")
  })

  it("works with verbose = FALSE without error", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")

    seu1 <- new_test_seurat(5L, 5L)
    seu2 <- new_test_seurat(5L, 5L, seed = 33L)

    suppressWarnings(
      merged <- SCIntegrate(seu1, seu2, method = "skip", verbose = FALSE)
    )

    expect_s4_class(merged, "Seurat")
  })
})

# get_names_4_ids ----------------------------------------------------------

describe("get_names_4_ids", {
  it("returns named arguments as-is", {
    result <- get_names_4_ids(a = 1, b = 2, c = 3)

    expect_equal(result, c("a", "b", "c"))
  })

  it("infers names from variable names for unnamed arguments", {
    foo <- 1
    bar <- 2

    result <- get_names_4_ids(foo, bar)

    expect_equal(result, c("foo", "bar"))
  })

  it("handles mixed named and unnamed arguments", {
    x <- 1
    result <- get_names_4_ids(a = 1, x, c = 3)

    expect_equal(result, c("a", "x", "c"))
  })

  it("returns character(0) for empty input", {
    result <- get_names_4_ids()

    expect_equal(result, character(0))
  })

  it("accepts .quoses parameter for unnamed args", {
    x <- 1
    y <- 2
    quos <- rlang::quos(x, y)
    result <- get_names_4_ids(x, y, .quoses = quos)

    expect_equal(result, c("x", "y"))
  })

  it("uses .quoses for unnamed arguments when length matches dots", {
    x <- 1
    # .quoses must have same length as dots (2 elements: a named, x unnamed)
    quos <- rlang::quos(a, x)
    result <- get_names_4_ids(a = 1, x, .quoses = quos)

    # Named arg "a" takes precedence; unnamed uses quosure label "x"
    expect_equal(result, c("a", "x"))
  })

  it("infers names correctly for single unnamed argument", {
    my_var <- 42
    result <- get_names_4_ids(my_var)

    expect_equal(result, "my_var")
  })
})
