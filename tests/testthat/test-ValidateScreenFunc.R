test_that("`ValidateScreenFunc` works", {
  func <- \(x) {
    z <- x + 1L
    # Chinese
    b
    s(z)
    ~`x`
    list(scRNA = 1L, z = z)
  }

  expect_error(ValidateScreenFunc(func))

  func2 <- \(x) {
    z <- x + 1L
    list(scRNA = 1L, z = z)
  }

  ValidateScreenFunc(func2)
})
