test_that("get_var_value works", {
  f <- function(save_path = "./analysis") {
    save_path_new <- file.path(save_path, "res")
    return(save_path_new)
  }
  f_res <- get_var_value("save_path_new", f) # returns "./analysis/res"

  expect_equal(f_res, "./analysis/res")

  g <- function(a = 1, b = 2) {
    c <- a * 2 + b * 3
    d = c^2
    e <<- d - 1
    e
  }
  g_res<-get_var_value("e", g) # returns 63

  expect_equal(g_res, 63)

  h <- function(a = "A", ...) {
    a <- 1
    return(a)
    a <- 2
    a
  }
  h_res<-get_var_value("a", h) # returns 1

  expect_equal(h_res, 1)
})
