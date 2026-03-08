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
})
