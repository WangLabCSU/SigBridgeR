test_that("setThreads", {
  setThreads(1L)

  setThreads(4L, tf_config = list(xla = T, inter_op = 4L, intra_op = 4L))
})
