test_that("`ValidateScreenFunc` works", {
  func <- \(x) {
    z <- x + 1
    # 中文
    b
    s(z)
    ~`×`
    list(scRNA_data = 1, z = z)
  }

  ValidateScreenFunc(func)
})
