# SymbolConvert -----------------------------------------------------------
# Real Ensembl IDs (grch37/grch38):
#   ENSG00000000003 -> TSPAN6, ENSG00000000419 -> DPM1, ENSG00000000457 -> SCYL3

new_test_mat <- function(
  ids = c(
    "ENSG00000000003",
    "ENSG00000000419",
    "ENSG00000000457"
  )
) {
  n <- length(ids)
  matrix(
    seq_len(n * 2),
    nrow = n,
    ncol = 2,
    dimnames = list(ids, c("S1", "S2"))
  )
}

describe("SymbolConvert - conversion", {
  it("converts matrix rownames to gene symbols", {
    mat <- new_test_mat()

    result <- SymbolConvert(mat, update_symbol = FALSE, verbose = FALSE)

    expect_identical(rownames(result), c("TSPAN6", "DPM1", "SCYL3"))
    expect_identical(colnames(result), c("S1", "S2"))
    expect_equal(dim(result), c(3, 2))
  })

  it("converts a vector of IDs to gene symbols", {
    result <- SymbolConvert(
      c("ENSG00000000003", "ENSG00000000419"),
      update_symbol = FALSE,
      verbose = FALSE
    )

    expect_identical(result, c("TSPAN6", "DPM1"))
  })

  it("strips version suffixes before conversion", {
    mat <- new_test_mat(c("ENSG00000000003.16", "ENSG00000000419.12"))

    result <- SymbolConvert(mat, update_symbol = FALSE, verbose = FALSE)

    expect_identical(rownames(result), c("TSPAN6", "DPM1"))
  })

  it("supports hg19 genome build", {
    result <- SymbolConvert(
      c("ENSG00000000003", "ENSG00000000419"),
      genome_build = "hg19",
      update_symbol = FALSE,
      verbose = FALSE
    )

    expect_identical(result, c("TSPAN6", "DPM1"))
  })
})

describe("SymbolConvert - NA handling", {
  it("replaces NA symbols with original IDs (version kept)", {
    mat <- new_test_mat(c("ENSG00000000003.16", "NOT_A_GENE"))

    expect_message(
      result <- SymbolConvert(mat, update_symbol = FALSE, verbose = FALSE),
      "NA values"
    )
    expect_identical(rownames(result), c("TSPAN6", "NOT_A_GENE"))
  })

  it("replaces NA symbols in vector input with original IDs", {
    expect_message(
      result <- SymbolConvert(
        c("ENSG00000000003", "NOT_A_GENE"),
        update_symbol = FALSE,
        verbose = FALSE
      ),
      "NA values"
    )
    expect_identical(result, c("TSPAN6", "NOT_A_GENE"))
  })

  it("does not message when no NA produced", {
    mat <- new_test_mat()

    expect_no_message(
      SymbolConvert(mat, update_symbol = FALSE, verbose = FALSE)
    )
  })
})

describe("SymbolConvert - duplicated symbols", {
  it("restores original IDs for duplicated symbols", {
    mat <- new_test_mat(c(
      "ENSG00000000003",
      "ENSG00000000003",
      "ENSG00000000419"
    ))

    expect_message(
      result <- SymbolConvert(mat, update_symbol = FALSE, verbose = FALSE),
      "duplicated gene symbols"
    )
    expect_identical(rownames(result), c("TSPAN6", "ENSG00000000003", "DPM1"))
  })

  it("does not message when no duplicates", {
    mat <- new_test_mat()

    expect_no_message(
      SymbolConvert(mat, update_symbol = FALSE, verbose = FALSE)
    )
  })
})

describe("SymbolConvert - update_symbol_method", {
  it("rejects invalid methods", {
    mat <- new_test_mat("ENSG00000000003")

    expect_error(
      SymbolConvert(mat, update_symbol_method = "bad"),
      "scCustomize"
    )
  })

  it("uses Seurat when method = 'Seurat'", {
    mat <- new_test_mat("ENSG00000000003")
    called <- new.env()
    testthat::local_mocked_bindings(
      UpdateSymbolList = function(x, verbose) {
        called$seurat <- TRUE
        paste0(x, "_seurat")
      },
      .package = "Seurat"
    )

    result <- SymbolConvert(
      mat,
      update_symbol_method = "Seurat",
      verbose = FALSE
    )

    expect_true(called$seurat)
    expect_identical(rownames(result), "TSPAN6_seurat")
  })

  it("uses scCustomize HGNC when method = 'scCustomize' and hg build", {
    mat <- new_test_mat("ENSG00000000003")
    called <- new.env()
    testthat::local_mocked_bindings(
      Updated_HGNC_Symbols = function(x, case_check_as_warn, verbose) {
        called$hgnc <- TRUE
        list(Output_Features = paste0(x, "_hgnc"))
      },
      .package = "scCustomize"
    )

    result <- SymbolConvert(
      mat,
      update_symbol_method = "scCustomize",
      verbose = FALSE
    )

    expect_true(called$hgnc)
    expect_identical(rownames(result), "TSPAN6_hgnc")
  })

  it("uses scCustomize MGI for mouse genome builds", {
    mat <- new_test_mat("ENSMUSG00000000001")
    testthat::local_mocked_bindings(
      convert_hm_genes = function(x, type, genome_build) "Trp53",
      .package = "IDConverter"
    )
    called <- new.env()
    testthat::local_mocked_bindings(
      Updated_MGI_Symbols = function(x, verbose) {
        called$mgi <- TRUE
        list(Output_Features = paste0(x, "_mgi"))
      },
      .package = "scCustomize"
    )

    result <- SymbolConvert(
      mat,
      genome_build = "mm10",
      update_symbol_method = "scCustomize",
      verbose = FALSE
    )

    expect_true(called$mgi)
    expect_identical(rownames(result), "Trp53_mgi")
  })

  it("skips symbol update when update_symbol = FALSE", {
    mat <- new_test_mat("ENSG00000000003")
    called <- new.env()
    testthat::local_mocked_bindings(
      UpdateSymbolList = function(x, verbose) {
        called$seurat <- TRUE
        x
      },
      .package = "Seurat"
    )

    result <- SymbolConvert(
      mat,
      update_symbol = FALSE,
      update_symbol_method = "Seurat",
      verbose = FALSE
    )

    expect_null(called$seurat)
    expect_identical(rownames(result), "TSPAN6")
  })
})

describe("SymbolConvert - deprecated unknown_format", {
  it("warns when unknown_format is supplied and ignores it", {
    mat <- new_test_mat(c("ENSG00000000003", "NOT_A_GENE"))

    expect_warning(
      expect_message(
        result <- SymbolConvert(
          mat,
          unknown_format = "unknown_{k}",
          update_symbol = FALSE,
          verbose = FALSE
        ),
        "NA values"
      ),
      class = "lifecycle_warning_deprecated"
    )
    expect_identical(rownames(result), c("TSPAN6", "NOT_A_GENE"))
  })

  it("does not warn when unknown_format is not supplied", {
    mat <- new_test_mat("ENSG00000000003")

    expect_no_warning(
      SymbolConvert(mat, update_symbol = FALSE, verbose = FALSE)
    )
  })
})

describe("SymbolConvert - input validation", {
  it("aborts on non-matrix-like, non-vector input", {
    expect_error(
      SymbolConvert(new.env(), update_symbol = FALSE, verbose = FALSE),
      "Ensembles version IDs or TCGA version IDs"
    )
  })

  it("aborts when matrix has no rownames", {
    mat <- matrix(1:4, nrow = 2)
    expect_error(
      SymbolConvert(mat, update_symbol = FALSE, verbose = FALSE),
      "Genes are missing"
    )
  })

  it("rejects invalid genome builds", {
    mat <- new_test_mat("ENSG00000000003")

    expect_error(
      SymbolConvert(mat, genome_build = "rn6", verbose = FALSE),
      "hg38"
    )
  })

  it("rejects non-flag update_symbol and verbose", {
    mat <- new_test_mat("ENSG00000000003")

    expect_error(
      SymbolConvert(mat, update_symbol = "yes", verbose = FALSE)
    )
    expect_error(
      SymbolConvert(mat, update_symbol = FALSE, verbose = 1L)
    )
  })
})
