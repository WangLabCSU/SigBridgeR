test_that("`ValidateScreenFunc` works", {
  func <- \(x) {
    z <- x + 1
    # Chinese
    b
    s(z)
    ~`x`
    list(scRNA = 1, z = z)
  }

  expect_error(ValidateScreenFunc(func))

  func2 <- \(x) {
    z <- x + 1
    list(scRNA = 1, z = z)
  }

  ValidateScreenFunc(func2)
})
