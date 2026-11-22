skip_if_not_installed("PhenoMapR")

# ---- local fixtures -------------------------------------------------------

.pm_bulk <- function(n_genes = 6L, n_samples = 4L) {
  matrix(
    seq_len(n_genes * n_samples),
    nrow = n_genes,
    dimnames = list(paste0("G", seq_len(n_genes)), paste0("S", seq_len(n_samples)))
  )
}

.pm_binary_pheno <- function(n = 4L) {
  setNames(rep(c(0, 1), length.out = n), paste0("S", seq_len(n)))
}

.pm_seurat <- function(n_cells = 4L) {
  counts <- Matrix::Matrix(
    matrix(seq_len(5L * n_cells), nrow = 5L, ncol = n_cells),
    sparse = TRUE
  )
  colnames(counts) <- paste0("C", seq_len(n_cells))
  suppressWarnings(SeuratObject::CreateSeuratObject(
    counts,
    min.cells = 0L,
    min.features = 0L
  ))
}

.pm_validate_args <- function(...) {
  modifyList(
    list(
      matched_bulk = .pm_bulk(),
      phenotype = .pm_binary_pheno(),
      label_type = "PhenoMapR",
      phenotype_class = "binary",
      verbose = FALSE,
      seed = 42L
    ),
    list(...)
  )
}

# ---- ValidatePhenoMapRParams: binary / continuous -------------------------

describe("ValidatePhenoMapRParams - binary and continuous phenotypes", {
  it("resolves defaults and normalizes a named phenotype vector", {
    pheno <- .pm_binary_pheno()
    p <- do.call(ValidatePhenoMapRParams, .pm_validate_args())

    expect_equal(p$label_type, "PhenoMapR")
    expect_equal(p$phenotype_class, "binary")
    expect_equal(p$phenotype_column, "phenotype_value")
    expect_null(p$survival_time)
    expect_null(p$survival_event)
    expect_equal(p$pheno_df$cell_id, names(pheno))
    expect_equal(p$pheno_df$phenotype_value, unname(pheno))
    expect_equal(p$normalize, TRUE)
    expect_equal(p$platform, c("auto", "rnaseq", "microarray"))
    expect_null(p$probe_annotation)
    expect_equal(p$hugo_species, c("human", "mouse"))
    expect_equal(p$binary_positive_reference, c("second", "first"))
    expect_equal(p$assay, "RNA")
    expect_equal(p$layer, "data")
    expect_equal(p$z_score_cutoff, 2)
    expect_equal(p$reference_sign, 1L)
    expect_equal(p$score_mode, c("weighted_sum", "activity_adjusted"))
    expect_equal(p$vars_to_regress, c("S.Score", "G2M.Score", "nCount_RNA"))
    expect_null(p$cancer_type)
    expect_false(p$pseudobulk)
    expect_null(p$group_by)
    expect_null(p$sample_metadata)
    expect_equal(p$percentile, 0.05)
    expect_null(p$permutation_n)
    expect_equal(p$permutation_seed, 42L)
    expect_false(p$verbose)
    expect_equal(p$seed, 42L)
    expect_equal(p$dots, list())
  })

  it("honours dots overrides for downstream PhenoMapR calls", {
    metadata <- data.frame(group = c("a", "b", "a", "b"))
    p <- do.call(
      ValidatePhenoMapRParams,
      .pm_validate_args(
        phenotype_column = "response",
        normalize = FALSE,
        platform = "microarray",
        probe_annotation = "annotation",
        hugo_species = "mouse",
        binary_positive_reference = "first",
        assay = "SCT",
        layer = "counts",
        z_score_cutoff = 3,
        reference_sign = -1L,
        score_mode = "activity_adjusted",
        vars_to_regress = "nFeature_RNA",
        cancer_type = "LUAD",
        pseudobulk = TRUE,
        group_by = "group",
        sample_metadata = metadata,
        percentile = 0.1,
        permutation_n = 50L,
        permutation_seed = 7L
      )
    )

    expect_equal(p$phenotype_column, "response")
    expect_equal(colnames(p$pheno_df), c("cell_id", "response"))
    expect_false(p$normalize)
    expect_equal(p$platform, "microarray")
    expect_equal(p$probe_annotation, "annotation")
    expect_equal(p$hugo_species, "mouse")
    expect_equal(p$binary_positive_reference, "first")
    expect_equal(p$assay, "SCT")
    expect_equal(p$layer, "counts")
    expect_equal(p$z_score_cutoff, 3)
    expect_equal(p$reference_sign, -1L)
    expect_equal(p$score_mode, "activity_adjusted")
    expect_equal(p$vars_to_regress, "nFeature_RNA")
    expect_equal(p$cancer_type, "LUAD")
    expect_true(p$pseudobulk)
    expect_equal(p$group_by, "group")
    expect_equal(p$sample_metadata, metadata)
    expect_equal(p$percentile, 0.1)
    expect_equal(p$permutation_n, 50L)
    expect_equal(p$permutation_seed, 7L)
    expect_equal(p$dots$normalize, FALSE)
  })

  it("uses the requested column name for continuous phenotypes", {
    pheno <- setNames(c(1.5, 2.5, 3.5, 4.5), paste0("S", 1:4))
    p <- do.call(
      ValidatePhenoMapRParams,
      .pm_validate_args(
        phenotype = pheno,
        phenotype_class = "continuous",
        phenotype_column = "dose"
      )
    )

    expect_equal(p$phenotype_class, "continuous")
    expect_equal(colnames(p$pheno_df), c("cell_id", "dose"))
    expect_equal(p$pheno_df$dose, unname(pheno))
  })

  it("builds a ScreenMethodConfig carrying the resolved parameters", {
    p <- do.call(ValidatePhenoMapRParams, .pm_validate_args(label_type = "my_label"))

    expect_true(S7::S7_inherits(p$cache_config, ScreenMethodConfig))
    expect_equal(p$cache_config@method_name, "PhenoMapR")
    expect_equal(p$cache_config@phenotype_class, "binary")
    expect_equal(p$cache_config@label_type, "my_label")
    expect_equal(p$cache_config@param$label_type, "my_label")
    expect_equal(p$cache_config@param$assay, "RNA")
  })
})

# ---- ValidatePhenoMapRParams: survival ------------------------------------

describe("ValidatePhenoMapRParams - survival phenotypes", {
  it("normalizes time/status columns into a cell_id/time/status data.frame", {
    pheno <- data.frame(
      time = c(1, 2, 3, 4),
      status = c(0, 1, 0, 1),
      row.names = paste0("S", 1:4)
    )
    p <- do.call(
      ValidatePhenoMapRParams,
      .pm_validate_args(
        phenotype = pheno,
        phenotype_class = "survival"
      )
    )

    expect_equal(p$phenotype_column, NULL)
    expect_equal(p$survival_time, "time")
    expect_equal(p$survival_event, "status")
    expect_equal(colnames(p$pheno_df), c("cell_id", "time", "status"))
    expect_equal(p$pheno_df$cell_id, rownames(pheno))
    expect_equal(p$pheno_df$time, pheno$time)
    expect_equal(p$pheno_df$status, pheno$status)
  })

  it("falls back to the first two columns for nonstandard column names", {
    pheno <- data.frame(
      followup = c(5, 6, 7, 8),
      death = c(1, 0, 1, 0),
      row.names = paste0("S", 1:4)
    )
    p <- do.call(
      ValidatePhenoMapRParams,
      .pm_validate_args(
        phenotype = pheno,
        phenotype_class = "survival"
      )
    )

    expect_equal(p$survival_time, "followup")
    expect_equal(p$survival_event, "death")
    expect_equal(p$pheno_df$time, pheno$followup)
    expect_equal(p$pheno_df$status, pheno$death)
  })

  it("respects explicit survival_time and survival_event overrides", {
    pheno <- data.frame(
      followup = c(5, 6, 7, 8),
      death = c(1, 0, 1, 0),
      other = c(9, 9, 9, 9),
      row.names = paste0("S", 1:4)
    )
    p <- do.call(
      ValidatePhenoMapRParams,
      .pm_validate_args(
        phenotype = pheno,
        phenotype_class = "survival",
        survival_time = "followup",
        survival_event = "death"
      )
    )

    expect_equal(p$survival_time, "followup")
    expect_equal(p$survival_event, "death")
    expect_equal(colnames(p$pheno_df), c("cell_id", "time", "status"))
  })
})

# ---- ValidatePhenoMapRParams: input validation ----------------------------

describe("ValidatePhenoMapRParams - input validation", {
  it("rejects a non-character or empty label_type", {
    expect_error(
      do.call(ValidatePhenoMapRParams, .pm_validate_args(label_type = 1L)),
      "label_type"
    )
    expect_error(
      do.call(ValidatePhenoMapRParams, .pm_validate_args(label_type = character(0L))),
      "length 1"
    )
  })

  it("rejects empty and unsupported phenotype_class values", {
    expect_error(
      do.call(ValidatePhenoMapRParams, .pm_validate_args(phenotype_class = character(0L))),
      "phenotype_class"
    )
    expect_error(
      do.call(ValidatePhenoMapRParams, .pm_validate_args(phenotype_class = "categorical")),
      "must be one of"
    )
  })

  it("requires sample IDs on the bulk matrix", {
    bulk <- .pm_bulk()
    colnames(bulk) <- NULL
    expect_error(
      do.call(ValidatePhenoMapRParams, .pm_validate_args(matched_bulk = bulk)),
      "sample IDs",
      fixed = TRUE
    )
  })

  it("requires names or row names on the phenotype", {
    expect_error(
      do.call(ValidatePhenoMapRParams, .pm_validate_args(phenotype = 0:1)),
      "must be named"
    )
  })

  it("requires phenotype IDs to match the bulk sample IDs", {
    pheno <- c(S1 = 0, S2 = 1, S3 = 0)
    expect_error(
      do.call(ValidatePhenoMapRParams, .pm_validate_args(phenotype = pheno)),
      "must match"
    )
  })

  it("requires a two-column data.frame for survival phenotypes", {
    expect_error(
      do.call(ValidatePhenoMapRParams, .pm_validate_args(phenotype_class = "survival")),
      "data.frame with time and status"
    )

    pheno <- data.frame(time = 1:4, row.names = paste0("S", 1:4))
    expect_error(
      do.call(
        ValidatePhenoMapRParams,
        .pm_validate_args(
          phenotype = pheno,
          phenotype_class = "survival"
        )
      ),
      "data.frame with time and status"
    )
  })
})

# ---- DoPhenoMapR ----------------------------------------------------------

describe("DoPhenoMapR", {
  it("maps PhenoMapR groups to the package-wide vocabulary", {
    bulk <- .pm_bulk()
    pheno <- .pm_binary_pheno()
    sc <- .pm_seurat()

    local_mocked_bindings(
      derive_reference_from_bulk = function(...) c(G1 = 1, G2 = -1),
      PhenoMap = function(...) {
        data.frame(
          cell_id = colnames(sc),
          score = seq_along(colnames(sc)),
          row.names = colnames(sc)
        )
      },
      define_phenotype_groups = function(...) {
        data.frame(
          cell_id = colnames(sc),
          group = c("Most Adverse", "Most Favorable", "Other", NA),
          row.names = colnames(sc)
        )
      },
      .package = "PhenoMapR"
    )

    res <- DoPhenoMapR(
      bulk,
      sc,
      pheno,
      phenotype_class = "binary",
      verbose = FALSE,
      seed = 42L
    )

    expect_named(res, c("scRNA_data", "phenomap_scores", "reference"))
    expect_s4_class(res$scRNA_data, "Seurat")
    expect_equal(res$reference, c(G1 = 1, G2 = -1))
    expect_equal(
      res$phenomap_scores$PhenoMapR,
      c("Positive", "Negative", "Neutral", NA_character_)
    )
    expect_equal(
      unname(res$scRNA_data$PhenoMapR),
      c("Positive", "Negative", "Neutral", NA_character_)
    )
    expect_true("score" %in% colnames(res$scRNA_data@meta.data))
    expect_equal(res$scRNA_data@misc$PhenoMapR$method_name, "PhenoMapR")
    expect_equal(res$scRNA_data@misc$PhenoMapR$label_type, "PhenoMapR")
  })

  it("forwards resolved parameters to the PhenoMapR functions", {
    bulk <- .pm_bulk()
    pheno <- .pm_binary_pheno()
    sc <- .pm_seurat()
    captured <- new.env(parent = emptyenv())

    local_mocked_bindings(
      derive_reference_from_bulk = function(...) {
        captured$derive <- list(...)
        c(G1 = 1)
      },
      PhenoMap = function(...) {
        captured$map <- list(...)
        data.frame(
          cell_id = colnames(sc),
          score = seq_along(colnames(sc)),
          row.names = colnames(sc)
        )
      },
      define_phenotype_groups = function(...) {
        captured$groups <- list(...)
        data.frame(
          cell_id = colnames(sc),
          group = rep("Other", ncol(sc)),
          row.names = colnames(sc)
        )
      },
      .package = "PhenoMapR"
    )

    DoPhenoMapR(
      bulk,
      sc,
      pheno,
      phenotype_class = "binary",
      label_type = "pm_label",
      normalize = FALSE,
      platform = "microarray",
      assay = "SCT",
      layer = "counts",
      z_score_cutoff = 3,
      reference_sign = -1L,
      score_mode = "activity_adjusted",
      vars_to_regress = "nFeature_RNA",
      percentile = 0.1,
      permutation_n = 50L,
      permutation_seed = 7L,
      verbose = FALSE,
      seed = 42L
    )

    expect_equal(captured$derive$bulk_expression, bulk)
    expect_equal(captured$derive$sample_id_column, "cell_id")
    expect_equal(captured$derive$phenotype_column, "phenotype_value")
    expect_equal(captured$derive$phenotype_type, "binary")
    expect_false(captured$derive$normalize)
    expect_equal(captured$derive$platform, "microarray")
    expect_true(inherits(captured$map$expression, "Seurat"))
    expect_equal(captured$map$reference, c(G1 = 1))
    expect_equal(captured$map$assay, "SCT")
    expect_equal(captured$map$layer, "counts")
    expect_equal(captured$map$z_score_cutoff, 3)
    expect_equal(captured$map$reference_sign, -1L)
    expect_equal(captured$map$score_mode, "activity_adjusted")
    expect_equal(captured$map$vars_to_regress, "nFeature_RNA")
    expect_equal(captured$map$permutation_n, 50L)
    expect_equal(captured$map$permutation_seed, 7L)
    expect_equal(captured$groups$percentile, 0.1)
    expect_null(captured$groups$score_columns)
  })

  it("passes survival time/status columns through to the reference builder", {
    bulk <- .pm_bulk()
    pheno <- nb_survival_phenotype(4L)
    sc <- .pm_seurat()
    captured <- new.env(parent = emptyenv())

    local_mocked_bindings(
      derive_reference_from_bulk = function(...) {
        captured$derive <- list(...)
        c(G1 = 1)
      },
      PhenoMap = function(...) {
        data.frame(
          cell_id = colnames(sc),
          score = seq_along(colnames(sc)),
          row.names = colnames(sc)
        )
      },
      define_phenotype_groups = function(...) {
        data.frame(
          cell_id = colnames(sc),
          group = rep("Other", ncol(sc)),
          row.names = colnames(sc)
        )
      },
      .package = "PhenoMapR"
    )

    DoPhenoMapR(
      bulk,
      sc,
      pheno,
      phenotype_class = "survival",
      verbose = FALSE,
      seed = 42L
    )

    expect_equal(captured$derive$phenotype_type, "survival")
    expect_equal(captured$derive$survival_time, "time")
    expect_equal(captured$derive$survival_event, "status")
    expect_null(captured$derive$phenotype_column)
    expect_equal(colnames(captured$derive$phenotype), c("cell_id", "time", "status"))
  })

  it("reports progress when verbose is TRUE", {
    bulk <- .pm_bulk()
    pheno <- .pm_binary_pheno()
    sc <- .pm_seurat()

    local_mocked_bindings(
      derive_reference_from_bulk = function(...) c(G1 = 1),
      PhenoMap = function(...) {
        data.frame(
          cell_id = colnames(sc),
          score = seq_along(colnames(sc)),
          row.names = colnames(sc)
        )
      },
      define_phenotype_groups = function(...) {
        data.frame(
          cell_id = colnames(sc),
          group = rep("Other", ncol(sc)),
          row.names = colnames(sc)
        )
      },
      .package = "PhenoMapR"
    )

    messages <- capture_messages(
      DoPhenoMapR(bulk, sc, pheno, phenotype_class = "binary", verbose = TRUE, seed = 42L)
    )
    all_messages <- paste(messages, collapse = "\n")

    expect_match(all_messages, "Start PhenoMapR screen", fixed = TRUE)
    expect_match(all_messages, "Deriving reference from bulk data", fixed = TRUE)
    expect_match(all_messages, "Scoring single cells", fixed = TRUE)
    expect_match(all_messages, "Defining phenotype groups", fixed = TRUE)
    expect_match(all_messages, "Done", fixed = TRUE)
  })

  it("propagates parameter validation errors", {
    expect_error(
      DoPhenoMapR(
        .pm_bulk(),
        .pm_seurat(),
        .pm_binary_pheno(),
        phenotype_class = "categorical",
        verbose = FALSE,
        seed = 42L
      ),
      "must be one of"
    )
  })
})
