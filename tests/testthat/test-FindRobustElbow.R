# Tests for FindRobustElbow() (R/11c-FindRobustElbow.R) and
# find_robust_elbow_cpp() (src/FindRobustElbow.cpp).
#
# Public contract:
#   * The per-PC statistics are obtained from the PCA cell embeddings
#     (cells x PCs), not from the `@stdev` slot, and are streamed through
#     tatami (`Rtatami.h`), so both base matrices and S4 `Matrix` classes are
#     accepted by `beachmat::initializeCpp()`.
#   * `ndims` is capped at the number of available PCs.
#   * The recommended dimension is the second largest value among the twelve
#     candidates of the six heuristics plus the fixed anchors 10/15/.../35.
#   * `FindRobustElbow()` aborts when the Seurat object has no `"pca"`
#     reduction, and validates `verbose` / `ndims` with chk.

# Reference implementation of the six heuristics, mirroring the original R code.
ref_find_robust_elbow <- function(embeddings, ndims = ncol(embeddings)) {
  embeddings <- as.matrix(embeddings)
  ndims <- min(ndims, ncol(embeddings))

  stdev <- apply(
    embeddings[, seq_len(ndims), drop = FALSE],
    2,
    stats::sd
  )
  variance <- stdev^2
  pct_variance <- variance / sum(variance) * 100
  cumulative_variance <- cumsum(pct_variance)

  first_above <- function(x, threshold, fallback) {
    idx <- which(x > threshold)[1L]
    if (is.na(idx)) fallback else idx
  }

  last_above <- function(x, threshold, fallback) {
    idx <- which(x > threshold)
    if (length(idx) > 0L) max(idx) else fallback
  }

  cumulative_90 <- first_above(cumulative_variance, 90, ndims)
  cumulative_80 <- first_above(cumulative_variance, 80, ndims)
  above_mean <- last_above(pct_variance, mean(pct_variance), 10L)
  above_2sd <- last_above(pct_variance, 2 * stats::sd(pct_variance), 10L)

  second_deriv <- if (ndims < 3L) {
    ndims
  } else {
    min(which.max(abs(diff(diff(variance)))) + 1L, ndims)
  }

  distance <- if (ndims < 2L) {
    1L
  } else {
    start <- c(1, variance[1L])
    end <- c(ndims, variance[ndims])
    distances <- vapply(
      seq_len(ndims),
      function(i) {
        point <- c(i, variance[i])
        abs(
          (end[2L] - start[2L]) *
            point[1L] -
            (end[1L] - start[1L]) * point[2L] +
            end[1L] * start[2L] -
            end[2L] * start[1L]
        ) /
          sqrt((end[2L] - start[2L])^2 + (end[1L] - start[1L])^2)
      },
      numeric(1L)
    )
    which.max(distances)
  }

  candidates <- c(
    cumulative_90,
    cumulative_80,
    above_mean,
    above_2sd,
    second_deriv,
    distance,
    10L,
    15L,
    20L,
    25L,
    30L,
    35L
  )

  list(
    stdev = stdev,
    pct_variance = pct_variance,
    cumulative_variance = cumulative_variance,
    cumulative_90 = cumulative_90,
    cumulative_80 = cumulative_80,
    above_mean = above_mean,
    above_2sd = above_2sd,
    second_deriv = second_deriv,
    distance = distance,
    # the fixed anchors may exceed the PCs that actually exist
    final_dims = min(sort(candidates, decreasing = TRUE)[2L], ndims)
  )
}

make_embeddings <- function(n_cells = 120L, n_pcs = 20L, seed = 42L) {
  withr::with_seed(
    seed,
    matrix(
      stats::rnorm(n_cells * n_pcs, sd = seq(5, 0.5, length.out = n_pcs)),
      nrow = n_cells,
      ncol = n_pcs
    )
  )
}

# find_robust_elbow_cpp -------------------------------------------------------

describe("find_robust_elbow_cpp", {
  it("reproduces the reference heuristics on a dense matrix", {
    embeddings <- make_embeddings()
    ref <- ref_find_robust_elbow(embeddings)
    res <- find_robust_elbow_cpp(beachmat::initializeCpp(embeddings), 20L)

    expect_equal(res$stdev, ref$stdev, tolerance = 1e-12)
    expect_equal(res$pct_variance, ref$pct_variance, tolerance = 1e-10)
    expect_equal(
      res$cumulative_variance,
      ref$cumulative_variance,
      tolerance = 1e-10
    )
    expect_equal(res$cumulative_90, ref$cumulative_90)
    expect_equal(res$cumulative_80, ref$cumulative_80)
    expect_equal(res$above_mean, ref$above_mean)
    expect_equal(res$above_2sd, ref$above_2sd)
    expect_equal(res$second_deriv, ref$second_deriv)
    expect_equal(res$distance, ref$distance)
    expect_equal(res$final_dims, ref$final_dims)
  })

  it("accepts S4 Matrix input with the same result", {
    embeddings <- make_embeddings()
    sparse <- methods::as(Matrix::Matrix(embeddings), "CsparseMatrix")

    dense_res <- find_robust_elbow_cpp(
      beachmat::initializeCpp(embeddings),
      20L
    )
    sparse_res <- find_robust_elbow_cpp(
      beachmat::initializeCpp(sparse),
      20L
    )

    expect_equal(sparse_res, dense_res, tolerance = 1e-12)
  })

  it("only uses the first `ndims` PCs and caps `ndims` at ncol", {
    embeddings <- make_embeddings()
    ref <- ref_find_robust_elbow(embeddings, ndims = 5L)
    res <- find_robust_elbow_cpp(beachmat::initializeCpp(embeddings), 5L)

    expect_length(res$pct_variance, 5L)
    expect_equal(res$final_dims, ref$final_dims)
    # the fixed anchors (10, 15, ...) cannot exceed the PCs that exist
    expect_lte(res$final_dims, 5L)

    capped <- find_robust_elbow_cpp(
      beachmat::initializeCpp(embeddings),
      1000L
    )
    expect_length(capped$pct_variance, ncol(embeddings))
    expect_lte(capped$final_dims, ncol(embeddings))
  })

  it("rejects inputs that cannot produce a dimension estimate", {
    expect_error(
      find_robust_elbow_cpp(beachmat::initializeCpp(make_embeddings()), 0L),
      "positive integer"
    )
    expect_error(
      find_robust_elbow_cpp(
        beachmat::initializeCpp(matrix(0, nrow = 1L, ncol = 5L)),
        5L
      ),
      "At least 2 cells"
    )
    expect_error(
      find_robust_elbow_cpp(
        beachmat::initializeCpp(matrix(0, nrow = 10L, ncol = 5L)),
        5L
      ),
      "zero variance"
    )

    non_finite <- make_embeddings()
    non_finite[1L, 1L] <- NaN
    expect_error(
      find_robust_elbow_cpp(beachmat::initializeCpp(non_finite), 20L),
      "non-finite"
    )
  })
})

# FindRobustElbow -------------------------------------------------------------

make_pca_seurat <- function(npcs = 10L) {
  withr::with_seed(1L, {
    counts <- matrix(
      stats::rpois(80L * 60L, 5),
      nrow = 80L,
      ncol = 60L
    )
    rownames(counts) <- paste0("gene", seq_len(80L))
    colnames(counts) <- paste0("cell", seq_len(60L))
  })

  obj <- suppressWarnings(SeuratObject::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "CsparseMatrix")
  ))
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = 50L, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  suppressWarnings(Seurat::RunPCA(obj, npcs = npcs, verbose = FALSE))
}

describe("FindRobustElbow", {
  it("matches the reference computed from the PCA stdev", {
    skip_if_not_installed("Seurat")
    skip_if_not_installed("SeuratObject")

    obj <- make_pca_seurat(npcs = 10L)
    # Seurat stores sd(embeddings) in the stdev slot, so both must agree
    ref <- ref_find_robust_elbow(obj[["pca"]]@cell.embeddings, ndims = 10L)

    expect_equal(
      FindRobustElbow(obj, verbose = FALSE, ndims = 10L),
      unname(ref$final_dims)
    )
  })

  it("is silent when verbose = FALSE and caps ndims", {
    skip_if_not_installed("Seurat")
    skip_if_not_installed("SeuratObject")

    obj <- make_pca_seurat(npcs = 10L)

    expect_no_message(
      res <- FindRobustElbow(obj, verbose = FALSE, ndims = 50L)
    )
    expect_type(res, "integer")
    expect_length(res, 1L)
    # anchors such as 15/20/... must not be recommended beyond the 10 PCs
    expect_lte(res, 10L)
  })

  it("reports the method results when verbose = TRUE", {
    skip_if_not_installed("Seurat")
    skip_if_not_installed("SeuratObject")

    obj <- make_pca_seurat(npcs = 10L)

    expect_message(
      FindRobustElbow(obj, verbose = TRUE, ndims = 10L),
      "Final Recommended Dimensions"
    )
    expect_warning(
      FindRobustElbow(obj, verbose = TRUE, ndims = 50L),
      "Only 10 PCs available"
    )
  })

  it("aborts when PCA has not been computed", {
    skip_if_not_installed("SeuratObject")

    obj <- suppressWarnings(SeuratObject::CreateSeuratObject(
      counts = matrix(1:20, nrow = 5L, ncol = 4L)
    ))

    expect_error(
      FindRobustElbow(obj, verbose = FALSE),
      "PCA has not been computed"
    )
  })

  it("validates verbose and ndims", {
    skip_if_not_installed("Seurat")
    skip_if_not_installed("SeuratObject")

    obj <- make_pca_seurat(npcs = 10L)

    expect_error(FindRobustElbow(obj, verbose = "yes"), class = "chk_error")
    expect_error(FindRobustElbow(obj, ndims = 0L), class = "chk_error")
  })
})

# plot_robust_elbow ----------------------------------------------------------

describe("plot_robust_elbow", {
  it("returns a ggplot with one layer per method candidate", {
    embeddings <- make_embeddings()
    res <- find_robust_elbow_cpp(beachmat::initializeCpp(embeddings), 20L)

    plot <- plot_robust_elbow(res)

    expect_s3_class(plot, "ggplot")
    # 1 point + 1 line + 6 vertical lines + 6 annotations
    expect_length(plot$layers, 14L)
    expect_equal(plot$data$PC, seq_along(res$pct_variance))
    expect_equal(plot$data$Variance, res$pct_variance)
  })
})
