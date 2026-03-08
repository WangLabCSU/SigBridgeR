test_that("It works", {
  res1 <- qs::qread("/home/data/sigbridger/scissor_result.qs", nthreads = 4L)
  res2 <- qs::qread("/home/data/sigbridger/scab_result.qs", nthreads = 4L)

  weight_test <- c("scissor" = 0.5, "scAB" = 0.5)

  merged <- MergeResult(res1, res2, weights = weight_test)
})
