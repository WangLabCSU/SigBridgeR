skip_if_not_installed("Seurat")

# Screen() - registry-driven dispatch --------------------------------------

test_that("GenerateRoxy.Screen.screen_method() lists every registered method", {
  lines <- GenerateRoxy.Screen.screen_method()

  expect_type(lines, "character")
  expect_match(lines[[1L]], "@param screen_method", fixed = TRUE)
  expect_setequal(
    grep("^ - ", lines, value = TRUE),
    paste0(" - \"", names(ScreenStrategy), "\": see also [Do", names(ScreenStrategy), "()]")
  )
})

test_that("Screen() rejects unknown methods and phenotype classes", {
  bulk <- nb_bulk(20L, 20L)
  sc <- nb_seurat(20L, 20L)
  pheno <- nb_binary_phenotype(20L)

  expect_error(Screen(bulk, sc, pheno, screen_method = "NoSuchMethod"))
  expect_error(
    Screen(
      bulk, sc, pheno,
      screen_method = "scAB",
      phenotype_class = "continuous"
    )
  )
})

test_that("Screen() applies the mapper, forwards dots and wraps the result", {
  bulk <- nb_bulk()
  sc <- nb_seurat()
  pheno <- nb_binary_phenotype()

  local({
    executor <- function(matched_bulk, sc_data, phenotype, label_type, phenotype_class, ...) {
      list(
        scRNA_data = sc_data,
        label = label_type,
        klass = phenotype_class,
        bulk_ncol = ncol(matched_bulk),
        pheno_len = length(phenotype),
        dots = rlang::list2(...)
      )
    }
    mapper <- function(params) {
      params$label_type <- paste0(params$label_type, "_mapped")
      params
    }

    ScreenStrategy$TestScreenMethod <- ScreenMethod(
      method_name = "TestScreenMethod",
      method_version = as.character(utils::packageVersion("Seurat")),
      executor = executor,
      phenotype_class = c("binary", "continuous"),
      mapper = mapper
    )
    withr::defer(rm("TestScreenMethod", envir = ScreenStrategy))

    res <- Screen(
      matched_bulk = bulk,
      sc_data = sc,
      phenotype = pheno,
      phenotype_class = "binary",
      screen_method = "TestScreenMethod",
      custom_arg = 42L
    )

    expect_true(S7::S7_inherits(res, ScreenMethodResult))
    expect_true(inherits(res@scRNA_data, "Seurat"))
    expect_equal(res@extra$label, "TestScreenMethod_mapped")
    expect_equal(res@extra$klass, "binary")
    expect_equal(res@extra$bulk_ncol, 100L)
    expect_equal(res@extra$pheno_len, 100L)
    expect_equal(res@extra$dots$custom_arg, 42L)
  })
})
