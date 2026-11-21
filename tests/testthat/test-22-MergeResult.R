test_that("MergeResult validates inputs", {
  expect_error(MergeResult(), "Input objects must be provided")
  expect_error(
    expect_warning(MergeResult(list(not_seurat = 1L), verbose = FALSE), "Skipping object"),
    "No valid Seurat objects"
  )
})

test_that("MergeResult combines metadata and keeps first duplicate column", {
  s1 <- nb_seurat(20L, 20L, seed = 1L)
  s2 <- nb_seurat(20L, 20L, seed = 2L)

  s1$Scissor <- rep(c("Positive", "Negative"), 10L)
  s2$scPP <- rep(c("Positive", "Neutral"), each = 10L)
  s1$shared <- "from_first"
  s2$shared <- "from_second"

  merged <- MergeResult(s1, s2, verbose = FALSE)

  expect_true(inherits(merged, "Seurat"))
  expect_equal(ncol(merged), 20L)
  expect_equal(unname(merged$shared), rep("from_first", 20L))
  expect_equal(unname(merged$Scissor), rep(c("Positive", "Negative"), 10L))
  expect_equal(unname(merged$scPP), rep(c("Positive", "Neutral"), each = 10L))
})

test_that("MergeResult accepts lists containing scRNA_data", {
  s1 <- nb_seurat(20L, 20L, seed = 1L)
  s2 <- nb_seurat(20L, 20L, seed = 2L)
  s2$scPAS <- rep("Positive", 20L)

  merged <- MergeResult(list(scRNA_data = s1, note = "ignored"), list(scRNA_data = s2))

  expect_equal(ncol(merged), 20L)
  expect_equal(unname(merged$scPAS), rep("Positive", 20L))
})

test_that("MergeResult performs an inner join on cell barcodes", {
  s1 <- nb_seurat(20L, 30L, seed = 1L)
  s2 <- nb_seurat(20L, 30L, seed = 2L)
  s2_small <- subset(s2, cells = paste0("C", 11:30))
  s2_small$scPP <- rep("Positive", 20L)

  expect_warning(
    merged <- MergeResult(s1, s2_small, verbose = FALSE),
    "intersection"
  )

  expect_equal(ncol(merged), 20L)
  expect_setequal(colnames(merged), paste0("C", 11:30))
  expect_equal(unname(merged$scPP), rep("Positive", 20L))
})

test_that("MergeResult adds a weighted vote column", {
  s1 <- nb_seurat(20L, 20L, seed = 1L)
  s2 <- nb_seurat(20L, 20L, seed = 2L)
  s1$Scissor <- c(rep("Positive", 10L), rep("Negative", 10L))
  s2$scPP <- c(rep("Positive", 10L), rep("Neutral", 10L))

  merged <- MergeResult(
    s1,
    s2,
    weights = c(Scissor = 2, scPP = 1),
    ties.method = "first",
    verbose = FALSE
  )

  expect_true("vote" %in% colnames(merged[[]]))
  expect_length(merged$vote, 20L)
  expect_setequal(unique(merged$vote), c("Positive", "Other"))
})

test_that("MergeResult merges misc entries from all inputs", {
  s1 <- nb_seurat(20L, 20L, seed = 1L)
  s2 <- nb_seurat(20L, 20L, seed = 2L)
  s1@misc$from_first <- 1L
  s2@misc$from_second <- 2L

  merged <- MergeResult(s1, s2, verbose = FALSE)

  expect_equal(merged@misc$from_first, 1L)
  expect_equal(merged@misc$from_second, 2L)
})

test_that("GetSlotNames reports available Seurat slots", {
  s <- nb_seurat(20L, 20L)

  expect_equal(GetSlotNames(s, "assays"), "RNA")
  expect_equal(GetSlotNames(s, "graphs"), NULL)
  expect_equal(GetSlotNames(s, "unknown"), NULL)
})

test_that("MergeSlot carries reductions from later objects", {
  s1 <- nb_seurat(20L, 30L, seed = 1L)
  s2 <- nb_seurat(20L, 30L, seed = 2L)
  embeddings <- matrix(
    stats::rnorm(60L),
    nrow = 30L,
    dimnames = list(colnames(s2), c("PC_1", "PC_2"))
  )
  s2[["pca"]] <- SeuratObject::CreateDimReducObject(embeddings = embeddings, key = "PC_")
  s2_small <- subset(s2, cells = paste0("C", 11:30))

  suppressWarnings(merged <- MergeResult(s1, s2_small, verbose = FALSE))

  expect_true("pca" %in% names(merged@reductions))
  expect_equal(rownames(merged@reductions$pca), paste0("C", 11:30))
})
