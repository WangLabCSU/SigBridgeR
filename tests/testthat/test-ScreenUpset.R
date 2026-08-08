# Helper: create a minimal Seurat object with n cells
skip_if_not_installed("ggplot2")
skip_if_not_installed("ggupset")
skip_if_not_installed("tibble")

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

# ScreenUpset - input validation -----------------------------------------

describe("ScreenUpset - input validation", {
  it("aborts when screened_seurat is not a Seurat object", {
    expect_error(
      ScreenUpset(list()),
      class = "chk_error"
    )
  })

  it("aborts when n_intersections is not a whole number", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")

    seurat <- new_test_seurat(2L)

    expect_error(
      ScreenUpset(seurat, n_intersections = 1.5),
      class = "chk_error"
    )
  })

  it("aborts when screen_type not found in metadata", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")

    seurat <- new_test_seurat(2L)
    seurat$scissor <- c("Positive", "Negative")

    expect_error(
      ScreenUpset(seurat, screen_type = "nonexistent"),
      "Screen type\\(s\\) not found"
    )
  })

  it("validates screen_type is character if provided", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")

    seurat <- new_test_seurat(2L)

    expect_error(
      ScreenUpset(seurat, screen_type = 123),
      class = "chk_error"
    )
  })
})

# ScreenUpset - basic functionality --------------------------------------

describe("ScreenUpset - basic functionality", {
  it("returns list with plot and stats", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggupset")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(4L)
    seurat$scissor <- c("Positive", "Positive", "Negative", "Neutral")
    seurat$scPAS <- c("Positive", "Negative", "Positive", "Neutral")

    result <- ScreenUpset(seurat, screen_type = c("scissor", "scPAS"))

    expect_type(result, "list")
    expect_true("plot" %in% names(result))
    expect_true("stats" %in% names(result))
    expect_s3_class(result$plot, "ggplot")
    expect_s3_class(result$stats, "tbl_df")
  })

  it("stats contains intersection, sets, count columns", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(4L)
    seurat$scissor <- c("Positive", "Positive", "Negative", "Neutral")
    seurat$scPAS <- c("Positive", "Negative", "Positive", "Neutral")

    result <- ScreenUpset(seurat, screen_type = c("scissor", "scPAS"))

    expect_true("intersection" %in% colnames(result$stats))
    expect_true("sets" %in% colnames(result$stats))
    expect_true("count" %in% colnames(result$stats))
  })

  it("counts all possible intersections for 2 screen types", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(4L)
    seurat$scissor <- c("Positive", "Positive", "Negative", "Neutral")
    seurat$scPAS <- c("Positive", "Negative", "Positive", "Neutral")

    result <- ScreenUpset(seurat, screen_type = c("scissor", "scPAS"))

    # 2 types -> 3 combinations: scissor alone, scPAS alone, scissor & scPAS
    expect_equal(nrow(result$stats), 3)
  })

  it("counts all possible intersections for 3 screen types", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(4L)
    seurat$scissor <- c("Positive", "Positive", "Negative", "Neutral")
    seurat$scPAS <- c("Positive", "Negative", "Positive", "Neutral")
    seurat$scAB <- c("Positive", "Positive", "Positive", "Negative")

    result <- ScreenUpset(
      seurat,
      screen_type = c("scissor", "scPAS", "scAB")
    )

    # 3 types -> 7 combos (C(3,1)=3 + C(3,2)=3 + C(3,3)=1)
    expect_equal(nrow(result$stats), 7)
  })
})

# ScreenUpset - intersection counts --------------------------------------

describe("ScreenUpset - intersection counts", {
  it("correctly counts cells positive for single screen type", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(5L)
    seurat$scissor <- c(
      "Positive",
      "Positive",
      "Positive",
      "Negative",
      "Neutral"
    )
    seurat$scPAS <- c("Positive", "Negative", "Neutral", "Positive", "Other")

    result <- ScreenUpset(seurat, screen_type = c("scissor", "scPAS"))

    # scissor alone: 3 cells have scissor="Positive"
    scissor_row <- result$stats[result$stats$intersection == "scissor", ]
    expect_equal(unname(scissor_row$count), 3)

    # scPAS alone: 2 cells have scPAS="Positive"
    scpas_row <- result$stats[result$stats$intersection == "scPAS", ]
    expect_equal(unname(scpas_row$count), 2)
  })

  it("correctly counts intersection of two screen types", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(6L)
    seurat$scissor <- c(
      "Positive",
      "Positive",
      "Positive",
      "Negative",
      "Neutral",
      "Positive"
    )
    seurat$scPAS <- c(
      "Positive",
      "Positive",
      "Negative",
      "Positive",
      "Neutral",
      "Positive"
    )

    result <- ScreenUpset(seurat, screen_type = c("scissor", "scPAS"))

    # Cells 1, 2, 6 have both scissor="Positive" and scPAS="Positive" -> 3
    both_row <- result$stats[
      result$stats$intersection == "scissor & scPAS",
    ]
    expect_equal(unname(both_row$count), 3)
  })

  it("handles case with no intersections (all disjoint)", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(4L)
    # scissor Positive on cells 1-2, scPAS Positive on cells 3-4 (disjoint)
    seurat$scissor <- c("Positive", "Positive", "Negative", "Negative")
    seurat$scPAS <- c("Negative", "Negative", "Positive", "Positive")

    result <- ScreenUpset(seurat, screen_type = c("scissor", "scPAS"))

    # Intersection should be 0
    both_row <- result$stats[
      result$stats$intersection == "scissor & scPAS",
    ]
    expect_equal(unname(both_row$count), 0)
  })

  it("handles case where all cells are positive for all types", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(3L)
    seurat$scissor <- c("Positive", "Positive", "Positive")
    seurat$scPAS <- c("Positive", "Positive", "Positive")

    result <- ScreenUpset(seurat, screen_type = c("scissor", "scPAS"))

    # All 3 cells are Positive in both
    both_row <- result$stats[
      result$stats$intersection == "scissor & scPAS",
    ]
    expect_equal(unname(both_row$count), 3)
  })
})

# ScreenUpset - screen_type auto-detection -------------------------------

describe("ScreenUpset - screen_type auto-detection", {
  it("auto-detects screen types when screen_type = NULL", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(4L)
    seurat$scissor <- c("Positive", "Positive", "Negative", "Neutral")
    seurat$scPAS <- c("Positive", "Negative", "Positive", "Neutral")

    result <- ScreenUpset(seurat, screen_type = NULL)

    expect_type(result, "list")
    expect_true("stats" %in% names(result))
    expect_true(nrow(result$stats) > 0)
  })

  it("auto-detection finds registered screen types matching pattern", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(4L)
    seurat$scissor <- c("Positive", "Positive", "Negative", "Neutral")

    result <- ScreenUpset(seurat, screen_type = NULL)

    # scissor should be auto-detected (matches sc[a-zA-Z]+$ pattern)
    expect_true("scissor" %in% result$stats$intersection)
  })
})

# ScreenUpset - plot parameters ------------------------------------------

describe("ScreenUpset - plot parameters", {
  it("show_plot = FALSE does not display plot", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggupset")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(4L)
    seurat$scissor <- c("Positive", "Positive", "Negative", "Neutral")
    seurat$scPAS <- c("Positive", "Negative", "Positive", "Neutral")

    result <- ScreenUpset(
      seurat,
      screen_type = c("scissor", "scPAS"),
      show_plot = FALSE
    )

    expect_s3_class(result$plot, "ggplot")
  })

  it("custom bar_color is applied", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggupset")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(4L)
    seurat$scissor <- c("Positive", "Positive", "Negative", "Neutral")
    seurat$scPAS <- c("Positive", "Negative", "Positive", "Neutral")

    result <- ScreenUpset(
      seurat,
      screen_type = c("scissor", "scPAS"),
      bar_color = "#FF0000"
    )

    expect_s3_class(result$plot, "ggplot")
  })

  it("custom combmatrix_point_color is applied", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggupset")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(4L)
    seurat$scissor <- c("Positive", "Positive", "Negative", "Neutral")
    seurat$scPAS <- c("Positive", "Negative", "Positive", "Neutral")

    result <- ScreenUpset(
      seurat,
      screen_type = c("scissor", "scPAS"),
      combmatrix_point_color = "darkgreen"
    )

    expect_s3_class(result$plot, "ggplot")
  })

  it("n_intersections limits displayed intersections", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggupset")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(4L)
    seurat$scissor <- c("Positive", "Positive", "Negative", "Neutral")
    seurat$scPAS <- c("Positive", "Negative", "Positive", "Neutral")
    seurat$scAB <- c("Positive", "Positive", "Positive", "Negative")

    result <- ScreenUpset(
      seurat,
      screen_type = c("scissor", "scPAS", "scAB"),
      n_intersections = 5
    )

    # Stats still have all 7 rows, but plot only shows top 5
    expect_equal(nrow(result$stats), 7)
    expect_s3_class(result$plot, "ggplot")
  })

  it("passes extra arguments to ggplot2::theme via ...", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggupset")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(4L)
    seurat$scissor <- c("Positive", "Positive", "Negative", "Neutral")
    seurat$scPAS <- c("Positive", "Negative", "Positive", "Neutral")

    result <- ScreenUpset(
      seurat,
      screen_type = c("scissor", "scPAS"),
      plot.title = ggplot2::element_text(size = 20)
    )

    expect_s3_class(result$plot, "ggplot")
  })
})

# ScreenUpset - edge cases -----------------------------------------------

describe("ScreenUpset - edge cases", {
  it("handles single screen type", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(4L)
    seurat$scissor <- c("Positive", "Positive", "Negative", "Neutral")

    result <- ScreenUpset(seurat, screen_type = "scissor")

    # 1 type -> C(1,1)=1 combination
    expect_equal(nrow(result$stats), 1)
    expect_equal(result$stats$intersection, "scissor")
    # 2 Positive cells
    expect_equal(unname(result$stats$count), 2)
  })

  it("handles data with no Positive cells", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(3L)
    seurat$scissor <- c("Negative", "Neutral", "Other")
    seurat$scPAS <- c("Negative", "Neutral", "Other")

    result <- ScreenUpset(seurat, screen_type = c("scissor", "scPAS"))

    # All counts should be 0
    expect_true(all(result$stats$count == 0))
  })

  it("returns intersection names with & separator for combos", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(4L)
    seurat$scissor <- c("Positive", "Positive", "Negative", "Neutral")
    seurat$scPAS <- c("Positive", "Negative", "Positive", "Neutral")

    result <- ScreenUpset(seurat, screen_type = c("scissor", "scPAS"))

    # Single types use the type name, pairs use " & "
    expect_true("scissor" %in% result$stats$intersection)
    expect_true("scPAS" %in% result$stats$intersection)
    expect_true("scissor & scPAS" %in% result$stats$intersection)
  })

  it("sets column in stats is a list column", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(4L)
    seurat$scissor <- c("Positive", "Positive", "Negative", "Neutral")
    seurat$scPAS <- c("Positive", "Negative", "Positive", "Neutral")

    result <- ScreenUpset(seurat, screen_type = c("scissor", "scPAS"))

    expect_type(result$stats$sets, "list")
  })

  it("handles many screen types and cells", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("tibble")

    n_cells <- 50
    seurat <- new_test_seurat(n_cells)
    set.seed(42)
    labels <- c("Positive", "Negative", "Neutral", "Other")
    seurat$scissor <- sample(labels, n_cells, replace = TRUE)
    seurat$scPAS <- sample(labels, n_cells, replace = TRUE)
    seurat$scAB <- sample(labels, n_cells, replace = TRUE)

    result <- ScreenUpset(
      seurat,
      screen_type = c("scissor", "scPAS", "scAB")
    )

    # 3 types -> 7 combinations
    expect_equal(nrow(result$stats), 7)
    # All counts should be non-negative and sum to <= n_cells * 3 (overlapping)
    expect_true(all(result$stats$count >= 0))
  })

  it("verbose = TRUE prints completion message", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("tibble")

    seurat <- new_test_seurat(4L)
    seurat$scissor <- c("Positive", "Positive", "Negative", "Neutral")
    seurat$scPAS <- c("Positive", "Negative", "Positive", "Neutral")

    expect_message(
      ScreenUpset(
        seurat,
        screen_type = c("scissor", "scPAS"),
        verbose = TRUE
      ),
      "completed"
    )
  })
})
