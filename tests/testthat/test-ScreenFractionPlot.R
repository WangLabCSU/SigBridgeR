# Helper: create a minimal Seurat object with n cells
skip_if_not_installed("ggplot2")
skip_if_not_installed("patchwork")

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

# complete_counts ---------------------------------------------------------

describe("complete_counts", {
  it("fills missing combinations with n = 0", {
    df <- data.frame(
      group = c("A", "A", "B"),
      status = c("Positive", "Negative", "Positive"),
      n = c(3, 2, 1)
    )
    result <- complete_counts(df, group, status)

    # Should have all 4 combos (A/B x Positive/Negative)
    expect_equal(nrow(result), 4)
    # B-Negative should be filled with n=0
    b_neg <- result[result$group == "B" & result$status == "Negative", ]
    expect_equal(b_neg$n, 0)
  })

  it("preserves existing counts", {
    df <- data.frame(
      group = c("A", "A"),
      status = c("Positive", "Negative"),
      n = c(5, 3)
    )
    result <- complete_counts(df, group, status)

    a_pos <- result[result$group == "A" & result$status == "Positive", ]
    expect_equal(a_pos$n, 5)
    a_neg <- result[result$group == "A" & result$status == "Negative", ]
    expect_equal(a_neg$n, 3)
  })

  it("handles single column completion", {
    df <- data.frame(
      group = c("A", "B"),
      n = c(1, 2)
    )
    result <- complete_counts(df, group)

    expect_equal(nrow(result), 2)
    expect_setequal(result$group, c("A", "B"))
  })

  it("handles three columns", {
    df <- data.frame(
      group = c("A", "A", "B"),
      status = c("Positive", "Positive", "Negative"),
      batch = c("X", "Y", "X"),
      n = c(1, 2, 3)
    )
    result <- complete_counts(df, group, status, batch)

    # complete() expands to all observed combos across all columns
    expect_true(nrow(result) >= nrow(df))
    # Existing values preserved
    a_px <- result[
      result$group == "A" & result$status == "Positive" & result$batch == "X",
    ]
    expect_equal(a_px$n, 1)
  })

  it("works without tidyr installed (expand.grid fallback)", {
    df <- data.frame(
      group = c("A", "B"),
      status = c("Positive", "Negative"),
      n = c(3, 1)
    )
    local_mocked_bindings(
      is_installed = function(...) FALSE,
      .package = "rlang"
    )
    result <- complete_counts(df, group, status)

    expect_equal(nrow(result), 4)
  })

  it("works with tidyr installed", {
    skip_if_not_installed("tidyr")
    df <- data.frame(
      group = c("A", "B"),
      status = c("Positive", "Negative"),
      n = c(3, 1)
    )
    result <- complete_counts(df, group, status)

    expect_equal(nrow(result), 4)
    a_neg <- result[result$group == "A" & result$status == "Negative", ]
    expect_equal(a_neg$n, 0)
  })

  it("handles empty input", {
    df <- data.frame(
      group = character(0),
      status = character(0),
      n = integer(0)
    )
    result <- complete_counts(df, group, status)

    expect_equal(nrow(result), 0)
  })
})

# ScreenFractionPlot - input validation ----------------------------------

describe("ScreenFractionPlot - input validation", {
  it("aborts when screened_seurat is not a Seurat object", {
    expect_error(
      ScreenFractionPlot(list(), group_by = "Source"),
      class = "chk_error"
    )
  })

  it("aborts when group_by is not a single character", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")

    seurat <- new_test_seurat(2L)

    expect_error(
      ScreenFractionPlot(seurat, group_by = 123),
      class = "chk_error"
    )
  })

  it("aborts when group_by has length > 1", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")

    seurat <- new_test_seurat(2L)

    expect_error(
      ScreenFractionPlot(seurat, group_by = c("a", "b")),
      class = "chk_error"
    )
  })

  it("aborts when group_by column not in metadata", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")

    seurat <- new_test_seurat(2L)

    expect_error(
      ScreenFractionPlot(seurat, group_by = "NonExistent"),
      "Grouping variable not found"
    )
  })

  it("aborts when screen_type not found in metadata", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")

    seurat <- new_test_seurat(2L)
    seurat$Source <- c("A", "B")

    expect_error(
      ScreenFractionPlot(seurat, screen_type = "nonexistent"),
      "Screen type\\(s\\) not found"
    )
  })

  it("validates plot_color is named if provided", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")

    seurat <- new_test_seurat(2L)

    expect_error(
      ScreenFractionPlot(seurat, plot_color = c("red", "blue")),
      class = "chk_error"
    )
  })
})

# ScreenFractionPlot - single screen type ---------------------------------

describe("ScreenFractionPlot - single screen type", {
  it("returns list with stats and plot", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")

    seurat <- new_test_seurat(3L)
    seurat$Source <- c("A", "A", "B")
    seurat$scissor <- c("Positive", "Negative", "Neutral")

    result <- ScreenFractionPlot(seurat, screen_type = "scissor")

    expect_type(result, "list")
    expect_true("stats" %in% names(result))
    expect_true("plot" %in% names(result))
    expect_s3_class(result$plot, "ggplot")
    expect_s3_class(result$stats, "data.frame")
  })

  it("stats data.frame contains expected columns", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")

    seurat <- new_test_seurat(3L)
    seurat$Source <- c("A", "A", "B")
    seurat$scissor <- c("Positive", "Negative", "Neutral")

    result <- ScreenFractionPlot(seurat, screen_type = "scissor")

    expect_true("Source" %in% colnames(result$stats))
    expect_true("scissor" %in% colnames(result$stats))
    expect_true("n" %in% colnames(result$stats))
    expect_true("Total" %in% colnames(result$stats))
    expect_true("Fraction" %in% colnames(result$stats))
  })

  it("calculates correct fractions", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")

    seurat <- new_test_seurat(4L)
    seurat$Source <- c("A", "A", "A", "B")
    seurat$scissor <- c("Positive", "Positive", "Negative", "Positive")

    result <- ScreenFractionPlot(seurat, screen_type = "scissor")

    stats <- result$stats
    # Group A: 2 Positive, 1 Negative -> 2/3, 1/3
    a_pos <- stats[stats$Source == "A" & stats$scissor == "Positive", ]
    expect_equal(a_pos$Fraction, 2 / 3, tolerance = 1e-8)
    a_neg <- stats[stats$Source == "A" & stats$scissor == "Negative", ]
    expect_equal(a_neg$Fraction, 1 / 3, tolerance = 1e-8)
    # Group B: 1 Positive -> 1
    b_pos <- stats[stats$Source == "B" & stats$scissor == "Positive", ]
    expect_equal(b_pos$Fraction, 1)
  })

  it("respects show_null = FALSE (default)", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")

    seurat <- new_test_seurat(4L)
    seurat$Source <- c("A", "A", "A", "B")
    # Group B only has Negative; no Positive in B
    seurat$scissor <- c("Positive", "Positive", "Negative", "Negative")

    result <- ScreenFractionPlot(
      seurat,
      screen_type = "scissor",
      show_null = FALSE
    )

    plot_data <- result$plot$data
    expect_true(all(plot_data$Fraction > 0))
  })

  it("show_null = TRUE includes zero-fraction groups", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")

    seurat <- new_test_seurat(4L)
    seurat$Source <- c("A", "A", "A", "B")
    seurat$scissor <- c("Positive", "Positive", "Negative", "Negative")

    result <- ScreenFractionPlot(
      seurat,
      screen_type = "scissor",
      show_null = TRUE
    )

    plot_data <- result$plot$data
    b_pos <- plot_data[
      plot_data$Source == "B" & plot_data$scissor == "Positive",
    ]
    expect_equal(nrow(b_pos), 1)
    expect_equal(b_pos$Fraction, 0)
  })

  it("show_plot = FALSE returns result without displaying", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")

    seurat <- new_test_seurat(2L)
    seurat$Source <- c("A", "B")
    seurat$scissor <- c("Positive", "Negative")

    result <- ScreenFractionPlot(
      seurat,
      screen_type = "scissor",
      show_plot = FALSE
    )

    expect_type(result, "list")
    expect_true("stats" %in% names(result))
  })

  it("uses custom plot_color", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")

    seurat <- new_test_seurat(2L)
    seurat$Source <- c("A", "B")
    seurat$scissor <- c("Positive", "Negative")

    custom_colors <- c(
      "Neutral" = "#AAAAAA",
      "Other" = "#AAAAAA",
      "Positive" = "#00FF00",
      "Negative" = "#0000FF"
    )
    result <- ScreenFractionPlot(
      seurat,
      screen_type = "scissor",
      plot_color = custom_colors
    )

    expect_s3_class(result$plot, "ggplot")
  })

  it("uses plot_title parameter", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")

    seurat <- new_test_seurat(2L)
    seurat$Source <- c("A", "B")
    seurat$scissor <- c("Positive", "Negative")

    result <- ScreenFractionPlot(
      seurat,
      screen_type = "scissor",
      plot_title = "My Custom Title"
    )

    expect_equal(result$plot$labels$title, "My Custom Title")
  })

  it("uses custom x_lab and y_lab", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")

    seurat <- new_test_seurat(2L)
    seurat$Source <- c("A", "B")
    seurat$scissor <- c("Positive", "Negative")

    result <- ScreenFractionPlot(
      seurat,
      screen_type = "scissor",
      x_lab = "My X",
      y_lab = "My Y"
    )

    expect_equal(result$plot$labels$x, "My X")
    expect_equal(result$plot$labels$y, "My Y")
  })
})

# ScreenFractionPlot - multiple screen types ------------------------------

describe("ScreenFractionPlot - multiple screen types", {
  it("returns stats list, plot list, and combined_plot", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("patchwork")

    seurat <- new_test_seurat(3L)
    seurat$Source <- c("A", "A", "B")
    seurat$scissor <- c("Positive", "Negative", "Neutral")
    seurat$scPAS <- c("Positive", "Positive", "Negative")

    result <- ScreenFractionPlot(
      seurat,
      screen_type = c("scissor", "scPAS")
    )

    expect_type(result, "list")
    expect_true("stats" %in% names(result))
    expect_true("plot" %in% names(result))
    expect_true("combined_plot" %in% names(result))
    expect_type(result$stats, "list")
    expect_equal(names(result$stats), c("scissor", "scPAS"))
    expect_type(result$plot, "list")
    expect_equal(names(result$plot), c("scissor", "scPAS"))
  })

  it("respects ncol for combined plot layout", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("patchwork")

    seurat <- new_test_seurat(3L)
    seurat$Source <- c("A", "A", "B")
    seurat$scissor <- c("Positive", "Negative", "Neutral")
    seurat$scPAS <- c("Positive", "Positive", "Negative")

    result <- ScreenFractionPlot(
      seurat,
      screen_type = c("scissor", "scPAS"),
      ncol = 1
    )

    expect_s3_class(result$combined_plot, "patchwork")
  })

  it("uses vectorized plot_title for multiple screen types", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")

    seurat <- new_test_seurat(3L)
    seurat$Source <- c("A", "A", "B")
    seurat$scissor <- c("Positive", "Negative", "Neutral")
    seurat$scPAS <- c("Positive", "Positive", "Negative")

    result <- ScreenFractionPlot(
      seurat,
      screen_type = c("scissor", "scPAS"),
      plot_title = c("Scissor Results", "scPAS Results")
    )

    expect_equal(result$plot$scissor$labels$title, "Scissor Results")
    expect_equal(result$plot$scPAS$labels$title, "scPAS Results")
  })

  it("auto-detects screen_type when screen_type = NULL", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")

    seurat <- new_test_seurat(3L)
    seurat$Source <- c("A", "A", "B")
    seurat$scissor <- c("Positive", "Negative", "Neutral")

    result <- ScreenFractionPlot(seurat, screen_type = NULL)

    expect_type(result, "list")
    expect_true("stats" %in% names(result))
  })
})

# ScreenFractionPlot - label_type from misc -------------------------------

describe("ScreenFractionPlot - label_type from misc", {
  it("uses label_type from Seurat misc slot if available", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")

    seurat <- new_test_seurat(2L)
    seurat$Source <- c("A", "B")
    seurat$scissor <- c("Positive", "Negative")
    seurat@misc$scissor_type <- "Scissor Label"

    result <- ScreenFractionPlot(seurat, screen_type = "scissor")

    expect_equal(result$plot$labels$fill, "Scissor Label")
  })

  it("falls back to screen_type name when misc slot missing", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")

    seurat <- new_test_seurat(2L)
    seurat$Source <- c("A", "B")
    seurat$scissor <- c("Positive", "Negative")

    result <- ScreenFractionPlot(seurat, screen_type = "scissor")

    expect_equal(result$plot$labels$fill, "scissor")
  })
})

# ScreenFractionPlot - edge cases ----------------------------------------

describe("ScreenFractionPlot - edge cases", {
  it("handles single group with single status", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")

    seurat <- new_test_seurat(2L)
    seurat$Source <- c("A", "A")
    seurat$scissor <- c("Positive", "Positive")

    result <- ScreenFractionPlot(seurat, screen_type = "scissor")

    stats <- result$stats
    a_pos <- stats[stats$Source == "A" & stats$scissor == "Positive", ]
    expect_equal(a_pos$Fraction, 1)
    expect_equal(a_pos$n, 2)
  })

  it("handles group with zero cells (all Neutral/Other)", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")

    seurat <- new_test_seurat(3L)
    seurat$Source <- c("A", "B", "B")
    seurat$scissor <- c("Neutral", "Positive", "Positive")

    result <- ScreenFractionPlot(seurat, screen_type = "scissor")

    stats <- result$stats
    # Group A: 1 Neutral -> Fraction for Positive=0, Negative=0
    a_pos <- stats[stats$Source == "A" & stats$scissor == "Positive", ]
    expect_equal(a_pos$Fraction, 0)
    # Group B: 2 Positive -> 1.0
    b_pos <- stats[stats$Source == "B" & stats$scissor == "Positive", ]
    expect_equal(b_pos$Fraction, 1)
  })

  it("passes extra arguments to ggplot2::theme via ...", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")

    seurat <- new_test_seurat(2L)
    seurat$Source <- c("A", "B")
    seurat$scissor <- c("Positive", "Negative")

    result <- ScreenFractionPlot(
      seurat,
      screen_type = "scissor",
      legend.position = "bottom"
    )

    expect_s3_class(result$plot, "ggplot")
  })
})
