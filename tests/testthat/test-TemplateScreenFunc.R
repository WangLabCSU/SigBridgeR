test_that("multiplication works", {
  TemplateScreenFunc(filename = "tests/testthat/my_screen.R")

  TemplateScreenFunc(
    filename = "tests/testthat/my_screen.R",
    documentation = TRUE
  )

    TemplateScreenFunc(
    filename = "current",
    documentation = TRUE
  )

      TemplateScreenFunc(
    filename = "current",
    open = FALSE
    documentation = TRUE
  )

    TemplateScreenFunc(
    filename = "tests/testthat/my_screen.R",
    append = FALSE,
    documentation = FALSE
  )
})
