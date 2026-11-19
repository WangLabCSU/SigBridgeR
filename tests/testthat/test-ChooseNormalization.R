# Tests for ChooseNormalization() and its helpers
# (R/11e-ChooseNormalization.R).
#
# Public contract:
#   * Named Seurat objects are compared after subsampling to `subset_size`.
#   * `metrics` reports raw diagnostics plus min-max normalized components,
#     a weighted `composite_score`, and a descending `rank`.
#   * `recommendation` is the method with the highest composite score.
#   * Marker retention is only computed when `known_hvgs` is supplied.
#   * Invalid method lists, cell counts, thresholds, and weights abort.

# Fixtures -----------------------------------------------------------------

# Two objects sharing the same counts but different normalized data:
#   * "good": normalized values have a mean independent of gene expression
#   * "bad": normalized values are log1p(counts), retaining mean-variance bias
make_normalization_obj <- function(
  n_genes = 150L,
  n_cells = 30L,
  seed = 1L,
  mode = c("good", "bad"),
  variable_features = NULL
) {
  mode <- match.arg(mode)
  withr::local_seed(seed)

  genes <- sprintf("g%03d", seq_len(n_genes))
  cells <- sprintf("c%02d", seq_len(n_cells))
  lambda <- rep(seq(0.5, 5, length.out = n_genes), times = n_cells)
  counts <- matrix(
    stats::rpois(n_genes * n_cells, lambda = lambda),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(genes, cells)
  )
  counts <- Matrix::Matrix(counts, sparse = TRUE)

  obj <- SeuratObject::CreateSeuratObject(
    counts = counts,
    min.cells = 0L,
    min.features = 0L
  )

  normalized <- if (mode == "good") {
    matrix(
      stats::rnorm(n_genes * n_cells, mean = 5, sd = 1),
      nrow = n_genes,
      ncol = n_cells,
      dimnames = list(genes, cells)
    )
  } else {
    log1p(as.matrix(counts))
  }

  SeuratObject::LayerData(obj, assay = "RNA", layer = "data") <- normalized
  Seurat::VariableFeatures(obj) <- variable_features %||% genes
  obj
}

marker_genes <- function() {
  list(
    T_cell = c("g001", "g002", "g003"),
    B_cell = c("g004", "g005")
  )
}

# ChooseNormalization ------------------------------------------------------

describe("ChooseNormalization - method comparison", {
  it("ranks methods and recommends the best composite score", {
    good <- make_normalization_obj(
      mode = "good",
      variable_features = c(
        "g001",
        "g002",
        "g003",
        "g006",
        "g007",
        "g008",
        "g009",
        "g010"
      )
    )
    bad <- make_normalization_obj(
      mode = "bad",
      variable_features = c(
        "g100",
        "g101",
        "g102",
        "g103",
        "g104",
        "g105",
        "g106",
        "g107"
      )
    )

    result <- ChooseNormalization(
      Good = good,
      Bad = bad,
      subset_size = 30L,
      known_hvgs = marker_genes(),
      verbose = FALSE
    )

    expect_named(result, c("metrics", "recommendation"))
    expect_s3_class(result$metrics, "data.table")
    expect_setequal(result$metrics$method, c("Good", "Bad"))
    expect_equal(result$metrics$rank, seq_len(nrow(result$metrics)))
    expect_equal(result$recommendation, result$metrics$method[1L])
    expect_true(all(diff(result$metrics$composite_score) <= 0))
    expect_equal(result$metrics$n_cells, rep(30L, 2L))
    expect_equal(result$metrics$n_genes, rep(150L, 2L))
    expect_equal(result$recommendation, "Good")

    # The good object retains known markers, the bad object does not.
    retention <- result$metrics$marker_retention[
      match(c("Good", "Bad"), result$metrics$method)
    ]
    expect_gt(retention[1L], retention[2L])
  })

  it("computes the composite score from normalized metrics and weights", {
    obj <- make_normalization_obj(mode = "good")
    weight <- c(
      variance_stability = 0.2,
      marker_signal = 0.5,
      dropout_robustness = 0.3
    )

    result <- ChooseNormalization(
      Only = obj,
      subset_size = 20L,
      known_hvgs = marker_genes(),
      weight = weight,
      verbose = FALSE
    )
    metrics <- result$metrics

    expected <- metrics$variance_stability *
      weight[["variance_stability"]] +
      metrics$marker_signal * weight[["marker_signal"]] +
      metrics$dropout_robustness * weight[["dropout_robustness"]]

    expect_equal(metrics$composite_score, expected)
    expect_equal(metrics$rank, 1L)
    expect_equal(result$recommendation, "Only")
  })

  it("returns neutral normalized scores for a single method", {
    obj <- make_normalization_obj(mode = "good")

    result <- ChooseNormalization(
      Only = obj,
      subset_size = 20L,
      known_hvgs = marker_genes(),
      verbose = FALSE
    )
    metrics <- result$metrics

    expect_false(is.na(metrics$marker_retention))
    expect_equal(metrics$variance_stability, 0.5)
    expect_equal(metrics$marker_signal, 0.5)
    expect_equal(metrics$dropout_robustness, 0.5)
    expect_equal(metrics$composite_score, 0.5)
  })

  it("keeps missing marker retention when known_hvgs is NULL", {
    obj <- make_normalization_obj(mode = "good")

    result <- ChooseNormalization(
      Only = obj,
      subset_size = 20L,
      known_hvgs = NULL,
      verbose = FALSE
    )
    metrics <- result$metrics

    expect_true(is.na(metrics$marker_retention))
    expect_true(is.na(metrics$marker_signal))
    expect_true(is.na(metrics$composite_score))
    expect_equal(metrics$rank, 1L)
    expect_equal(result$recommendation, "Only")
  })

  it("subsets cells and is reproducible for a fixed seed", {
    good <- make_normalization_obj(mode = "good")
    bad <- make_normalization_obj(mode = "bad")

    first <- ChooseNormalization(
      Good = good,
      Bad = bad,
      subset_size = 10L,
      seed = 42L,
      verbose = FALSE
    )
    second <- ChooseNormalization(
      Good = good,
      Bad = bad,
      subset_size = 10L,
      seed = 42L,
      verbose = FALSE
    )

    expect_equal(first$metrics$n_cells, rep(10L, 2L))
    expect_equal(first$metrics, second$metrics)
  })

  it("defaults subset_size to the number of available cells", {
    obj <- make_normalization_obj(mode = "good")

    skip_if_not_installed("patchwork")

    expect_message(
      result <- ChooseNormalization(
        Only = obj,
        subset_size = integer(),
        verbose = TRUE
      ),
      "Using 30 cells"
    )

    expect_equal(result$metrics$n_cells, 30L)
    expect_equal(result$recommendation, "Only")
  })

  it("falls back to the package verbose and seed options", {
    obj <- make_normalization_obj(mode = "good")

    withr::local_options(list(
      SigBridgeR.verbose = FALSE,
      SigBridgeR.seed = 123L
    ))

    result <- ChooseNormalization(
      Only = obj
    )

    expect_equal(result$metrics$n_cells, 30L)
    expect_equal(result$recommendation, "Only")
  })

  it("reports the top method and key differentiators when verbose", {
    good <- make_normalization_obj(
      mode = "good",
      variable_features = c("g001", "g002", "g003", "g006", "g007")
    )
    bad <- make_normalization_obj(
      mode = "bad",
      variable_features = c("g100", "g101", "g102", "g103", "g104")
    )

    skip_if_not_installed("patchwork")

    expect_message(
      result <- ChooseNormalization(
        Good = good,
        Bad = bad,
        subset_size = 30L,
        known_hvgs = marker_genes(),
        verbose = TRUE
      ),
      "Method Ranking"
    )
    expect_equal(result$recommendation, "Good")
  })
})

# ChooseNormalizationCheck -------------------------------------------------

describe("ChooseNormalizationCheck", {
  it("accepts a valid method list", {
    obj <- make_normalization_obj(n_cells = 10L)

    expect_null(
      ChooseNormalizationCheck(
        method_objects = list(A = obj),
        method_names = "A",
        n_cells = 10L
      )
    )
  })

  it("requires between one and thirty methods", {
    obj <- make_normalization_obj(n_cells = 10L)

    expect_error(
      ChooseNormalizationCheck(list(), character(), integer()),
      "length"
    )

    many <- rep(list(obj), 31L)
    names(many) <- paste0("M", seq_along(many))
    expect_error(
      ChooseNormalizationCheck(
        method_objects = many,
        method_names = names(many),
        n_cells = rep(10L, length(many))
      ),
      "length"
    )
  })

  it("requires named methods", {
    obj <- make_normalization_obj(n_cells = 10L)

    expect_error(
      ChooseNormalizationCheck(
        method_objects = list(obj),
        method_names = "",
        n_cells = 10L
      ),
      "Options must be named"
    )
  })

  it("requires Seurat objects", {
    expect_error(
      ChooseNormalizationCheck(
        method_objects = list(A = 1L),
        method_names = "A",
        n_cells = 1L
      ),
      "Seurat"
    )
  })

  it("requires normalized data in the data layer", {
    counts <- Matrix::Matrix(
      matrix(
        1L:20L,
        nrow = 4L,
        dimnames = list(
          paste0("g", 1L:4L),
          paste0("c", 1L:5L)
        )
      ),
      sparse = TRUE
    )
    # CreateAssay5Object stores only the provided counts layer.
    assay <- SeuratObject::CreateAssay5Object(counts = counts)
    obj <- SeuratObject::CreateSeuratObject(counts = assay)

    expect_error(
      ChooseNormalizationCheck(
        method_objects = list(A = obj),
        method_names = "A",
        n_cells = 5L
      ),
      "normalized data"
    )
    expect_error(
      ChooseNormalizationCheck(
        method_objects = list(A = obj),
        method_names = "A",
        n_cells = 5L,
        layer = "data"
      ),
      "normalized data"
    )
  })

  it("requires matching cell counts", {
    small <- make_normalization_obj(n_cells = 10L)
    large <- make_normalization_obj(n_cells = 11L)

    expect_error(
      ChooseNormalizationCheck(
        method_objects = list(A = small, B = large),
        method_names = c("A", "B"),
        n_cells = c(10L, 11L)
      ),
      "different cell counts"
    )
  })

  it("rejects an invalid low-expression threshold", {
    obj <- make_normalization_obj(n_cells = 10L)

    expect_error(
      ChooseNormalizationCheck(
        method_objects = list(A = obj),
        method_names = "A",
        n_cells = 10L,
        low_expressed_thresh = 1.5
      ),
      "range|between"
    )
  })

  it("validates the weight vector", {
    obj <- make_normalization_obj(n_cells = 10L)
    check <- function(weight) {
      ChooseNormalizationCheck(
        method_objects = list(A = obj),
        method_names = "A",
        n_cells = 10L,
        weight = weight
      )
    }

    expect_error(
      check(c(0.4, 0.35, 0.25)),
      "named|names"
    )
    expect_error(
      check(c(foo = 0.4, bar = 0.35, baz = 0.25)),
      "weight must contain exactly"
    )
    expect_error(
      check(c(variance_stability = 0.5, marker_signal = 0.5)),
      "length"
    )
    expect_error(
      check(c(
        variance_stability = 0.4,
        marker_signal = 0.35,
        dropout_robustness = 0.3
      )),
      "weight must sum to 1"
    )
  })
})

# ExtractMetrics -----------------------------------------------------------

describe("ExtractMetrics", {
  it("returns raw diagnostics without marker retention", {
    obj <- make_normalization_obj(
      n_cells = 20L,
      variable_features = paste0("g", sprintf("%03d", 1L:10L))
    )

    result <- ExtractMetrics(
      obj = obj,
      assay = "RNA",
      method_name = "Only",
      ground_truth_markers = NULL
    )

    expect_named(
      result,
      c(
        "method",
        "variance_mean_cor",
        "marker_retention",
        "mean_dropout_residual",
        "n_cells",
        "n_genes",
        "n_hvgs_detected"
      )
    )
    expect_equal(result$method, "Only")
    expect_true(is.finite(result$variance_mean_cor))
    expect_true(is.finite(result$mean_dropout_residual))
    expect_true(is.na(result$marker_retention))
    expect_equal(result$n_cells, 20L)
    expect_equal(result$n_genes, 150L)
    expect_equal(result$n_hvgs_detected, 10L)
  })

  it("computes marker retention against the first n_hvgs variable features", {
    obj <- make_normalization_obj(
      n_cells = 20L,
      variable_features = c("g001", "g002", "g010", "g020")
    )

    result <- ExtractMetrics(
      obj = obj,
      assay = "RNA",
      method_name = "Only",
      ground_truth_markers = marker_genes(),
      n_hvgs = 4L
    )

    # g001 and g002 are markers; g010 and g020 are not.
    expect_equal(result$marker_retention, 0.5)
    expect_equal(result$n_hvgs_detected, 4L)
  })

  it("returns NA retention when no variable features are detected", {
    obj <- make_normalization_obj(
      n_cells = 20L,
      variable_features = character(0L)
    )

    result <- ExtractMetrics(
      obj = obj,
      assay = "RNA",
      method_name = "Only",
      ground_truth_markers = marker_genes()
    )

    expect_true(is.na(result$marker_retention))
    expect_equal(result$n_hvgs_detected, 0L)
  })

  it("returns NA retention when n_hvgs selects no genes", {
    obj <- make_normalization_obj(
      n_cells = 20L,
      variable_features = c("g001", "g002")
    )

    result <- ExtractMetrics(
      obj = obj,
      assay = "RNA",
      method_name = "Only",
      ground_truth_markers = marker_genes(),
      n_hvgs = 0L
    )

    expect_true(is.na(result$marker_retention))
    expect_equal(result$n_hvgs_detected, 2L)
  })
})

# ChooseNormalizationViz ---------------------------------------------------

describe("ChooseNormalizationViz", {
  it("builds a patchwork visualization from the metrics table", {
    skip_if_not_installed("patchwork")

    metrics <- data.table::data.table(
      method = c("A", "B"),
      variance_mean_cor = c(0.1, 0.5),
      composite_score = c(0.9, 0.4)
    )

    plot <- ChooseNormalizationViz(metrics)

    expect_s3_class(plot, "patchwork")
  })
})
