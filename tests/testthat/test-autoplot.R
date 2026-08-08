skip_if_not_installed("ggplot2")

# Helper: create a minimal Seurat object with n cells
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
  suppressWarnings(SeuratObject::CreateSeuratObject(counts = m))
}

# palette_SigBridgeR ------------------------------------------------------

describe("palette_SigBridgeR", {
  it("returns a character vector", {
    result <- palette_SigBridgeR()
    expect_type(result, "character")
    expect_true(length(result) > 0L)
  })

  it("all elements are valid hex color codes", {
    result <- palette_SigBridgeR()
    expect_match(result, "^#[0-9A-Fa-f]{6}$", all = TRUE)
  })

  it("returns exactly 40 colors", {
    result <- palette_SigBridgeR()
    expect_length(result, 40L)
  })

  it("returns all unique colors", {
    result <- palette_SigBridgeR()
    expect_equal(length(unique(tolower(result))), length(result))
  })
})

# autoplot.Seurat ---------------------------------------------------------

describe("autoplot.Seurat", {
  it("errors when required packages are not installed", {
    skip_if_not_installed("SeuratObject")

    seurat <- new_test_seurat()

    local_mocked_bindings(
      check_installed = function(
        pkg,
        reason = NULL,
        call = rlang::caller_env()
      ) {
        rlang::abort("Mock: package not installed", call = call)
      },
      .package = "SigBridgeR"
    )

    expect_error(
      autoplot(seurat),
      "Mock: package not installed"
    )
  })

  it("returns a ggplot object", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("SeuratObject")

    seurat <- new_test_seurat()

    local_mocked_bindings(
      check_installed = function(pkg, ...) invisible(),
      .package = "SigBridgeR"
    )
    local_mocked_bindings(
      DimPlot = function(object, reduction, ...) ggplot2::ggplot(),
      NoAxes = function(...) ggplot2::theme(),
      .package = "Seurat"
    )
    local_mocked_bindings(
      theme_dr = function(...) ggplot2::theme(),
      .package = "tidydr"
    )

    p <- autoplot(seurat)
    expect_s3_class(p, "ggplot")
  })

  it("uses custom cols when provided", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("SeuratObject")

    seurat <- new_test_seurat()
    custom_cols <- c("#FF0000", "#00FF00", "#0000FF")

    captured <- new.env(parent = emptyenv())

    local_mocked_bindings(
      check_installed = function(pkg, ...) invisible(),
      .package = "SigBridgeR"
    )
    local_mocked_bindings(
      DimPlot = function(object, reduction, cols, ...) {
        captured$cols <- cols
        ggplot2::ggplot()
      },
      NoAxes = function(...) ggplot2::theme(),
      .package = "Seurat"
    )
    local_mocked_bindings(
      theme_dr = function(...) ggplot2::theme(),
      .package = "tidydr"
    )

    autoplot(seurat, cols = custom_cols)
    expect_equal(captured$cols, custom_cols)
  })

  it("uses palette_SigBridgeR as default colors when cols is NULL", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("SeuratObject")

    seurat <- new_test_seurat()

    captured <- new.env(parent = emptyenv())

    local_mocked_bindings(
      check_installed = function(pkg, ...) invisible(),
      .package = "SigBridgeR"
    )
    local_mocked_bindings(
      DimPlot = function(object, reduction, cols, ...) {
        captured$cols <- cols
        ggplot2::ggplot()
      },
      NoAxes = function(...) ggplot2::theme(),
      .package = "Seurat"
    )
    local_mocked_bindings(
      theme_dr = function(...) ggplot2::theme(),
      .package = "tidydr"
    )

    autoplot(seurat)
    expect_equal(captured$cols, palette_SigBridgeR())
  })

  it("passes group.by to Seurat::DimPlot", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("SeuratObject")

    seurat <- new_test_seurat()

    captured <- new.env(parent = emptyenv())

    local_mocked_bindings(
      check_installed = function(pkg, ...) invisible(),
      .package = "SigBridgeR"
    )
    local_mocked_bindings(
      DimPlot = function(object, reduction, group.by, ...) {
        captured$group.by <- group.by
        ggplot2::ggplot()
      },
      NoAxes = function(...) ggplot2::theme(),
      .package = "Seurat"
    )
    local_mocked_bindings(
      theme_dr = function(...) ggplot2::theme(),
      .package = "tidydr"
    )

    autoplot(seurat, group.by = "orig.ident")
    expect_equal(captured$group.by, "orig.ident")
  })

  it("passes label and label.size to Seurat::DimPlot", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("SeuratObject")

    seurat <- new_test_seurat()

    captured <- new.env(parent = emptyenv())

    local_mocked_bindings(
      check_installed = function(pkg, ...) invisible(),
      .package = "SigBridgeR"
    )
    local_mocked_bindings(
      DimPlot = function(object, reduction, label, label.size, ...) {
        captured$label <- label
        captured$label.size <- label.size
        ggplot2::ggplot()
      },
      NoAxes = function(...) ggplot2::theme(),
      .package = "Seurat"
    )
    local_mocked_bindings(
      theme_dr = function(...) ggplot2::theme(),
      .package = "tidydr"
    )

    autoplot(seurat, label = FALSE, label.size = 6L)
    expect_false(captured$label)
    expect_equal(captured$label.size, 6L)
  })

  it("passes pt.size to Seurat::DimPlot", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("SeuratObject")

    seurat <- new_test_seurat()

    captured <- new.env(parent = emptyenv())

    local_mocked_bindings(
      check_installed = function(pkg, ...) invisible(),
      .package = "SigBridgeR"
    )
    local_mocked_bindings(
      DimPlot = function(object, reduction, pt.size, ...) {
        captured$pt.size <- pt.size
        ggplot2::ggplot()
      },
      NoAxes = function(...) ggplot2::theme(),
      .package = "Seurat"
    )
    local_mocked_bindings(
      theme_dr = function(...) ggplot2::theme(),
      .package = "tidydr"
    )

    autoplot(seurat, pt.size = 1.5)
    expect_equal(captured$pt.size, 1.5)
  })

  it("forwards ... to Seurat::DimPlot", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("SeuratObject")

    seurat <- new_test_seurat()

    captured <- new.env(parent = emptyenv())

    local_mocked_bindings(
      check_installed = function(pkg, ...) invisible(),
      .package = "SigBridgeR"
    )
    local_mocked_bindings(
      DimPlot = function(object, reduction, ...) {
        dots <- list(...)
        captured$dots <- dots
        ggplot2::ggplot()
      },
      NoAxes = function(...) ggplot2::theme(),
      .package = "Seurat"
    )
    local_mocked_bindings(
      theme_dr = function(...) ggplot2::theme(),
      .package = "tidydr"
    )

    autoplot(seurat, raster = TRUE, order = TRUE)
    expect_true(captured$dots$raster)
    expect_true(captured$dots$order)
  })

  it("uses default label and label.size when not specified", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("SeuratObject")

    seurat <- new_test_seurat()

    captured <- new.env(parent = emptyenv())

    local_mocked_bindings(
      check_installed = function(pkg, ...) invisible(),
      .package = "SigBridgeR"
    )
    local_mocked_bindings(
      DimPlot = function(object, reduction, label, label.size, ...) {
        captured$label <- label
        captured$label.size <- label.size
        ggplot2::ggplot()
      },
      NoAxes = function(...) ggplot2::theme(),
      .package = "Seurat"
    )
    local_mocked_bindings(
      theme_dr = function(...) ggplot2::theme(),
      .package = "tidydr"
    )

    autoplot(seurat)
    expect_false(captured$label)
    expect_equal(captured$label.size, 4)
  })

  it("uses 'umap' as default reduction", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("SeuratObject")

    seurat <- new_test_seurat()

    captured <- new.env(parent = emptyenv())

    local_mocked_bindings(
      check_installed = function(pkg, ...) invisible(),
      .package = "SigBridgeR"
    )
    local_mocked_bindings(
      DimPlot = function(object, reduction, ...) {
        captured$reduction <- reduction
        ggplot2::ggplot()
      },
      NoAxes = function(...) ggplot2::theme(),
      .package = "Seurat"
    )
    local_mocked_bindings(
      theme_dr = function(...) ggplot2::theme(),
      .package = "tidydr"
    )

    autoplot(seurat)
    expect_equal(captured$reduction, "umap")
  })

  it("passes reduction to Seurat::DimPlot", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("SeuratObject")

    seurat <- new_test_seurat()

    captured <- new.env(parent = emptyenv())

    local_mocked_bindings(
      check_installed = function(pkg, ...) invisible(),
      .package = "SigBridgeR"
    )
    local_mocked_bindings(
      DimPlot = function(object, reduction, ...) {
        captured$reduction <- reduction
        ggplot2::ggplot()
      },
      NoAxes = function(...) ggplot2::theme(),
      .package = "Seurat"
    )
    local_mocked_bindings(
      theme_dr = function(...) ggplot2::theme(),
      .package = "tidydr"
    )

    autoplot(seurat, reduction = "pca")
    expect_equal(captured$reduction, "pca")
  })

  it("uses ScreenStrategy colors when group.by is a registered strategy", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("SeuratObject")

    seurat <- new_test_seurat()

    captured <- new.env(parent = emptyenv())

    local_mocked_bindings(
      check_installed = function(pkg, ...) invisible(),
      .package = "SigBridgeR"
    )
    local_mocked_bindings(
      DimPlot = function(object, reduction, cols, ...) {
        captured$cols <- cols
        ggplot2::ggplot()
      },
      NoAxes = function(...) ggplot2::theme(),
      .package = "Seurat"
    )
    local_mocked_bindings(
      theme_dr = function(...) ggplot2::theme(),
      .package = "tidydr"
    )

    autoplot(seurat, group.by = "Scissor")
    expect_equal(
      captured$cols,
      c(
        "Other" = "#CECECE",
        "Neutral" = "#CECECE",
        "Positive" = "#c24b4b",
        "Negative" = "#5189bb"
      )
    )
  })

  it("uses palette_SigBridgeR when group.by is not a registered strategy", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("SeuratObject")

    seurat <- new_test_seurat()

    captured <- new.env(parent = emptyenv())

    local_mocked_bindings(
      check_installed = function(pkg, ...) invisible(),
      .package = "SigBridgeR"
    )
    local_mocked_bindings(
      DimPlot = function(object, reduction, cols, ...) {
        captured$cols <- cols
        ggplot2::ggplot()
      },
      NoAxes = function(...) ggplot2::theme(),
      .package = "Seurat"
    )
    local_mocked_bindings(
      theme_dr = function(...) ggplot2::theme(),
      .package = "tidydr"
    )

    autoplot(seurat, group.by = "orig.ident")
    expect_equal(captured$cols, palette_SigBridgeR())
  })
})

# autoplot.SigBridgeR::ScreenMethodResult ---------------------------------

describe("autoplot.SigBridgeR::ScreenMethodResult", {
  it("returns a ggplot object", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("SeuratObject")

    seurat <- new_test_seurat()
    screen_result <- ScreenMethodResult(scRNA_data = seurat)

    local_mocked_bindings(
      check_installed = function(pkg, ...) invisible(),
      .package = "SigBridgeR"
    )
    local_mocked_bindings(
      DimPlot = function(object, reduction, ...) ggplot2::ggplot(),
      NoAxes = function(...) ggplot2::theme(),
      .package = "Seurat"
    )
    local_mocked_bindings(
      theme_dr = function(...) ggplot2::theme(),
      .package = "tidydr"
    )

    p <- autoplot(screen_result)
    expect_s3_class(p, "ggplot")
  })

  it("delegates to autoplot.Seurat with scRNA_data", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("SeuratObject")

    seurat <- new_test_seurat()
    screen_result <- ScreenMethodResult(scRNA_data = seurat)

    captured <- new.env(parent = emptyenv())

    local_mocked_bindings(
      check_installed = function(pkg, ...) invisible(),
      .package = "SigBridgeR"
    )
    local_mocked_bindings(
      DimPlot = function(object, reduction, ...) {
        captured$seurat_object <- object
        ggplot2::ggplot()
      },
      NoAxes = function(...) ggplot2::theme(),
      .package = "Seurat"
    )
    local_mocked_bindings(
      theme_dr = function(...) ggplot2::theme(),
      .package = "tidydr"
    )

    autoplot(screen_result)
    expect_identical(captured$seurat_object, seurat)
  })

  it("passes extra arguments through to autoplot.Seurat", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("SeuratObject")

    seurat <- new_test_seurat()
    screen_result <- ScreenMethodResult(scRNA_data = seurat)

    captured <- new.env(parent = emptyenv())

    local_mocked_bindings(
      check_installed = function(pkg, ...) invisible(),
      .package = "SigBridgeR"
    )
    local_mocked_bindings(
      DimPlot = function(object, reduction, cols, label, ...) {
        captured$cols <- cols
        captured$label <- label
        ggplot2::ggplot()
      },
      NoAxes = function(...) ggplot2::theme(),
      .package = "Seurat"
    )
    local_mocked_bindings(
      theme_dr = function(...) ggplot2::theme(),
      .package = "tidydr"
    )

    custom_cols <- c("#FF0000", "#0000FF")
    autoplot(screen_result, cols = custom_cols, label = FALSE)
    expect_equal(captured$cols, custom_cols)
    expect_false(captured$label)
  })

  it("validates scRNA_data is a Seurat object at construction", {
    skip_if_not_installed("SeuratObject")

    expect_error(
      ScreenMethodResult(scRNA_data = "not_a_seurat"),
      "must be a.*Seurat"
    )
  })
})
