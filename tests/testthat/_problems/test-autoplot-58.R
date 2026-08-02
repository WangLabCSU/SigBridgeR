# Extracted from test-autoplot.R:58

# prequel ----------------------------------------------------------------------
new_test_seurat <- function(n_cells = 2L, n_genes = 2L) {
  m <- matrix(
    seq_len(n_cells * n_genes),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(
      paste0("gene", seq_len(n_genes)),
      paste0("cell", seq_len(n_cells))
    )
  )
  SeuratObject::CreateSeuratObject(counts = m)
}

# test -------------------------------------------------------------------------
it("errors when required packages are not installed", {
    skip_if_not_installed("SeuratObject")

    seurat <- new_test_seurat()

    local_mocked_bindings(
      check_installed = function(pkg, reason = NULL, call = rlang::caller_env()) {
        rlang::abort("Mock: package not installed", call = call)
      },
      .package = "rlang"
    )

    expect_error(
      autoplot(seurat),
      "Mock: package not installed"
    )
  })
