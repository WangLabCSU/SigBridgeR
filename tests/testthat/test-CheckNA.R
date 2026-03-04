test_that("it works", {
  # Example 1: Vector with NAs
  vec <- c(10, 20, NA, 40, NA, 60)
  result <- CheckNA(vec)
  # Output: Found 2 NA values
  #         Positions: 3, 5

  # Example 2: Named vector
  named_vec <- c(Sample1 = 10, Sample2 = NA, Sample3 = 30, Sample4 = NA)
  result <- CheckNA(named_vec)
  # Output: Found 2 NA values
  #         Positions: 2, 4
  #         with names: Sample2, Sample4

  # Example 3: Data frame with row and column names
  df <- data.frame(
    Gene1 = c(1.0, 2.0, NA),
    Gene2 = c(NA, 2.5, 3.0),
    Gene3 = c(1.5, NA, 3.5),
    row.names = c("Cell1", "Cell2", "Cell3")
  )
  result <- CheckNA(df)
  # Output: Found 3 NA values in data
  #         First 3 positions:
  #         Row 3 (Cell3), Col 1 (Gene1)
  #         Row 1 (Cell1), Col 2 (Gene2)
  #         Row 2 (Cell2), Col 3 (Gene3)

  # Example 4: Matrix without names
  mat <- matrix(c(1, NA, 3, 4, 5, NA), nrow = 2, ncol = 3)
  result <- CheckNA(mat)
  # Output: Found 2 NA values in data
  #         First 2 positions:
  #         Row 2, Col 1
  #         Row 1, Col 3
  expect_equal(result$count, 2)
  expect_equal(dim(result$positions), c(2, 2)) 

  # Example 5: No NAs found
  clean_data <- c(1, 2, 3, 4, 5)
  result <- CheckNA(clean_data)
  # Output: No NA values found in the vector
  expect_equal(result$count, 0)
  expect_equal(result$positions, NULL)

  # Example 6: Programmatic use - extract NA count
  df <- data.frame(a = c(1, NA, 3), b = c(NA, 2, 3))
  na_info <- CheckNA(df)
  expect_equal(na_info$count, 2) # Returns: 2
  na_info$positions # Returns: data frame with positions
})
