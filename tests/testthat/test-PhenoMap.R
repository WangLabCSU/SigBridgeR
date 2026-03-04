test_that("vector works", {
  v <- mtcars$mpg

  v2 <- PhenoMap(v, v > 10 ~ 1, v <= 10 ~ 0)

  expect_equal(v2, as.numeric(v > 10))

  a <- mtcars$cyl

  a2 <- PhenoMap(a, a > 4 ~ 1, a <= 4 ~ 0)

  expect_equal(a2, as.numeric(a > 4))

  f <- mtcars$hp

  f2 <- PhenoMap(f, f > 110 ~ 1, f < 110 ~ 0, f == 110 ~ 0.5)

  expect_equal(f2, as.numeric(ifelse(f > 110, 1, ifelse(f < 110, 0, 0.5))))

  f3 <- PhenoMap(f, f > 110 ~ "high", f < 110 ~ "low", f == 110 ~ "mid")

  expect_equal(f3, ifelse(f > 110, "high", ifelse(f < 110, "low", "mid")))
})

test_that("vector works", {
  b <- 1

  b2 <- PhenoMap(b, b > 1 ~ 1, b < 1 ~ 0)

  expect_equal(b2, NA)
})

test_that("data works", {
  c <- PhenoMap(mtcars, mpg > 15 ~ 1, mpg <= 15 ~ 0)

  expect_equal(mtcars$mpg, as.numeric(mtcars$mpg > 15))

  e <- PhenoMap(mtcars, mpg > 15 ~ "high", mpg <= 15 ~ "low", .default = -1)

  expect_equal(
    mtcars$mpg,
    ifelse(mtcars$mpg > 15, "high", ifelse(mtcars$mpg <= 15, "low", -1))
  )

  expect_error(PhenoMap(
    mtcars,
    mpg > 15 ~ "high",
    mpg <= 15 ~ 0,
    .default = TRUE
  ))
})

test_that(".default works", {
  c <- PhenoMap(mtcars, hp > 110 ~ 1, hp < 110 ~ 0, .default = -1)

  i1 <- which(c$hp == -1)
  i2 <- which(mtcars$hp == 110)

  expect_equal(i1, i2)
})
