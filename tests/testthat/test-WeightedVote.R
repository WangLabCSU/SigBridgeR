test_that("It works", {
  seurat <- qs::qread(
    "/home/data/sigbridger/merged_without_pipet.qs",
    nthreads = 4L
  )
  meta <- seurat[[]]

  vote_data <- meta[c("scAB", "scissor", "scPAS", "scPP", "DEGAS", "LP_SGL")]
  voter_weights <- setNames(
    abs(rnorm(6, mean = 0, sd = 0.5)),
    c("scAB", "scissor", "scPAS", "scPP", "DEGAS", "LP_SGL")
  )

  res <- WeightedVote(vote_data, voter_weights)

  expect_equal(length(res), nrow(vote_data))

  # compatibility test with data.table
  vote_data2 <- data.table::as.data.table(vote_data)
  res <- WeightedVote(vote_data2, voter_weights)
})
