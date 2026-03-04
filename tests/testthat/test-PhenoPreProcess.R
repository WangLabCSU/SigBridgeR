test_that("multiplication works", {
  # Example 1: Binary phenotype
  bulk_data <- matrix(rpois(100, 10), nrow = 20, ncol = 5)
  colnames(bulk_data) <- paste0("Sample", 1:5)
  rownames(bulk_data) <- paste0("Gene", 1:20)

  pheno_binary <- c(
    Sample1 = 1,
    Sample2 = 0,
    Sample3 = 1,
    Sample4 = 0,
    Sample5 = 1
  )

  result <- PhenoPreProcess(
    bulk = bulk_data,
    phenotype = pheno_binary,
    phenotype_class = "binary"
  )

  # Example 2: Continuous phenotype with discretization
  pheno_age <- c(
    Sample1 = 25,
    Sample2 = 35,
    Sample3 = 45,
    Sample4 = 55,
    Sample5 = 65
  )

  result <- PhenoPreProcess(
    bulk = bulk_data,
    phenotype = pheno_age,
    phenotype_class = "binary",
    pheno_age < 40 ~ 0,
    pheno_age >= 40 ~ 1
  )

  # Example 3: Survival data with automatic column detection
  pheno_surv <- data.frame(
    time = c(12, 24, 18, 36, 30),
    status = c(1, 0, 1, 1, 0),
    row.names = paste0("Sample", 1:5)
  )

  result <- PhenoPreProcess(
    bulk = bulk_data,
    phenotype = pheno_surv,
    phenotype_class = "survival"
  )

  # Example 4: Survival data with explicit column names
  pheno_surv_custom <- data.frame(
    follow_up = c(12, 24, 18, 36, 30),
    event = c(1, 0, 1, 1, 0),
    row.names = paste0("Sample", 1:5)
  )

  result <- PhenoPreProcess(
    bulk = bulk_data,
    phenotype = pheno_surv_custom,
    phenotype_class = "survival",
    select = c("follow_up", "event")
  )
})
