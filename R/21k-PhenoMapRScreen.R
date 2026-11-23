# ---- 2. DO PHENOMAPR ----

#' Validate and Prepare DoPhenoMapR Parameters
#'
#' @description
#' Internal helper that handles package checks, sample-id alignment, phenotype
#' normalization, input validation, and default-value resolution for
#' [DoPhenoMapR()].
#'
#' @inheritParams Screen
#' @param ... Additional dots forwarded from [DoPhenoMapR()].
#'
#' @return A named list with the normalized phenotype data frame (`pheno_df`),
#'   the derived parameters, and a `ScreenMethodConfig` (`cache_config`).
#'
#' @keywords internal
#' @family PhenoMapR
ValidatePhenoMapRParams <- function(
  matched_bulk,
  phenotype,
  label_type,
  phenotype_class,
  verbose,
  seed,
  ...
) {
  rlang::check_installed("PhenoMapR", action = \(pkg, ...) {
    rlang::check_installed("pak")
    pak::pak("brooksbenard/PhenoMapR")
  })

  chk::chk_character(label_type)
  chk::chk_length(label_type)
  chk::chk_length(phenotype_class)
  phenotype_class <- rlang::arg_match0(
    phenotype_class,
    c("binary", "survival", "continuous")
  )

  # -- sample-id alignment --------------------------------------------------
  sample_ids <- colnames(matched_bulk)
  if (is.null(sample_ids)) {
    rlang::abort("`matched_bulk` must have sample IDs in `colnames()`.")
  }

  pheno_ids <- if (is.data.frame(phenotype)) {
    rownames(phenotype)
  } else {
    names(phenotype)
  }
  if (is.null(pheno_ids)) {
    Abort(
      "`phenotype` must be named, or have row names for survival data frames."
    )
  }
  if (!setequal(pheno_ids, sample_ids)) {
    Abort(c(
      "Sample IDs of `phenotype` must match `colnames(matched_bulk)`.",
      "Keep only bulk samples that are linked to a phenotype."
    ))
  }

  dots <- rlang::list2(...)

  # -- normalize phenotype to the data.frame expected by PhenoMapR ----------
  if (phenotype_class == "survival") {
    if (!is.data.frame(phenotype) || ncol(phenotype) < 2L) {
      Abort(
        "`phenotype` must be a data.frame with time and status columns when `phenotype_class = \"survival\"`."
      )
    }
    survival_time <- dots$survival_time %||% "time"
    survival_event <- dots$survival_event %||% "status"
    if (!survival_time %in% colnames(phenotype)) {
      survival_time <- colnames(phenotype)[[1L]]
    }
    if (!survival_event %in% colnames(phenotype)) {
      survival_event <- colnames(phenotype)[[2L]]
    }
    pheno_df <- data.frame(
      cell_id = pheno_ids,
      time = phenotype[[survival_time]],
      status = phenotype[[survival_event]],
      stringsAsFactors = FALSE
    )
    phenotype_column <- NULL
  } else {
    phenotype_column <- dots$phenotype_column %||% "phenotype_value"
    pheno_df <- data.frame(
      cell_id = pheno_ids,
      value = unname(phenotype),
      stringsAsFactors = FALSE
    )
    names(pheno_df)[2L] <- phenotype_column
    survival_time <- NULL
    survival_event <- NULL
  }

  # -- defaults -------------------------------------------------------------
  res <- list(
    label_type = label_type,
    phenotype_class = phenotype_class,
    phenotype_column = phenotype_column,
    survival_time = survival_time,
    survival_event = survival_event,
    normalize = dots$normalize %||% TRUE,
    platform = dots$platform %||% c("auto", "rnaseq", "microarray"),
    probe_annotation = dots$probe_annotation,
    hugo_species = dots$hugo_species %||% c("human", "mouse"),
    binary_positive_reference = dots$binary_positive_reference %||% c("second", "first"),
    assay = dots$assay %||% "RNA",
    layer = dots$layer %||% "data",
    z_score_cutoff = dots$z_score_cutoff %||% 2L,
    reference_sign = dots$reference_sign %||% 1L,
    score_mode = dots$score_mode %||% c("weighted_sum", "activity_adjusted"),
    vars_to_regress = dots$vars_to_regress %||% c("S.Score", "G2M.Score", "nCount_RNA"),
    cancer_type = dots$cancer_type,
    pseudobulk = dots$pseudobulk %||% FALSE,
    group_by = dots$group_by,
    sample_metadata = dots$sample_metadata,
    percentile = dots$percentile %||% 0.05,
    permutation_n = dots$permutation_n,
    permutation_seed = dots$permutation_seed %||% seed,
    verbose = verbose,
    seed = seed,
    dots = dots
  )

  cache_config <- ScreenMethodConfig(
    method_name = "PhenoMapR",
    param = res,
    phenotype_class = phenotype_class,
    label_type = label_type
  )

  c(res, list(pheno_df = pheno_df, cache_config = cache_config))
}

#' @title PhenoMapR: Bulk-to-Single-Cell Phenotype Mapping
#'
#' @description
#' Derives a gene-level reference z-score from matched bulk expression and a
#' phenotype (binary, continuous, or survival), then scores single cells and
#' assigns tail phenotype groups via PhenoMapR.
#'
#' @param matched_bulk Matrix or data frame of bulk expression (genes x samples).
#'   Column names must match `phenotype` names/row names.
#' @param sc_data A Seurat object containing the cells to be scored.
#' @param phenotype Phenotype data. A named vector for binary/continuous traits,
#'   or a data frame with time-to-event and event columns for survival traits;
#'   names/row names must match `colnames(matched_bulk)`.
#' @param label_type Character. `meta.data` column receiving the phenotype group.
#'   Default: `"PhenoMapR"`.
#' @param phenotype_class One of `"binary"`, `"survival"`, or `"continuous"`.
#' @inheritParams PhenoMapR::derive_reference_from_bulk
#' @inheritParams PhenoMapR::PhenoMap
#' @inheritParams PhenoMapR::define_phenotype_groups
#' @param ... Additional arguments reserved for method compatibility.
#' @param seed Integer. Random seed. Default: `getFuncOption("seed")`.
#' @param verbose Logical. Whether to print progress messages.
#'   Default: `getFuncOption("verbose")`.
#'
#' @details
#' `define_phenotype_groups()` labels the upper tail of the score distribution
#' `"Most Adverse"` and the lower tail `"Most Favorable"`. The score increases
#' towards the phenotype-valued-as-1 end when the reference is derived with the
#' default orientation (for a numeric `0`/`1` phenotype,
#' `binary_positive_reference = "second"` regresses on the `"1"` level, and for
#' survival traits the Cox z-score is positive for hazard-increasing genes).
#' Therefore, when phenotype `1` denotes the unfavourable direction, the
#' `"Most Adverse"` tail is the one to be labelled `"Positive"`; the group
#' labels are mapped to the package-wide vocabulary as
#' `"Most Adverse" -> "Positive"`, `"Most Favorable" -> "Negative"` and
#' `"Other" -> "Neutral"`. Setting `reference_sign = -1` mirrors the score axis,
#' in which case the two tails swap. If a phenotype is encoded so that the
#' adverse class is *not* the reference-positive level (e.g.
#' `binary_positive_reference = "first"` with numeric `0`/`1` labels), the
#' score direction is reversed as well; choose the arguments so that higher
#' scores track the adverse phenotype.
#'
#' @return A list containing:
#'   - `scRNA_data`: Seurat object with `label_type` (`"Positive"`,
#'     `"Negative"`, `"Neutral"`, or `NA`) and a `PhenoMapR_*` score column in
#'     `meta.data`.
#'   - `phenomap_scores`: Per-cell scores from PhenoMap.
#'   - `reference`: Gene-level reference z-scores.
#'
#' @references
#' PhenoMapR: mapping phenotypic information from bulk to single-cell data.
#' <https://github.com/brooksbenard/PhenoMapR>
#'
#' @examples
#' \dontrun{
#' res <- DoPhenoMapR(
#'   matched_bulk = bulk_matrix,
#'   sc_data = seurat_obj,
#'   phenotype = a_named_vector,
#'   label_type = "PhenoMapR",
#'   phenotype_class = "binary"
#' )
#' }
#'
#' @export
#' @family screen_method
#' @family PhenoMapR
DoPhenoMapR <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = "PhenoMapR",
  phenotype_class = c("binary", "survival", "continuous"),
  phenotype_column = NULL,
  survival_time = "time",
  survival_event = "status",
  normalize = TRUE,
  platform = c("auto", "rnaseq", "microarray"),
  probe_annotation = NULL,
  hugo_species = c("human", "mouse"),
  binary_positive_reference = c("second", "first"),
  assay = "RNA",
  layer = "data",
  z_score_cutoff = 2L,
  reference_sign = 1L,
  score_mode = c("weighted_sum", "activity_adjusted"),
  vars_to_regress = c("S.Score", "G2M.Score", "nCount_RNA"),
  cancer_type = NULL,
  pseudobulk = FALSE,
  group_by = NULL,
  sample_metadata = NULL,
  percentile = 0.05,
  permutation_n = NULL,
  permutation_seed = NULL,
  ...,
  seed = SigBridgeRUtils::getFuncOption("seed"),
  verbose = SigBridgeRUtils::getFuncOption("verbose")
) {
  set.seed(seed)

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green("Start PhenoMapR screen"))
  }

  p <- ValidatePhenoMapRParams(
    matched_bulk = matched_bulk,
    phenotype = phenotype,
    label_type = label_type,
    phenotype_class = phenotype_class,
    phenotype_column = phenotype_column,
    survival_time = survival_time,
    survival_event = survival_event,
    normalize = normalize,
    platform = platform,
    probe_annotation = probe_annotation,
    hugo_species = hugo_species,
    binary_positive_reference = binary_positive_reference,
    assay = assay,
    layer = layer,
    z_score_cutoff = z_score_cutoff,
    reference_sign = reference_sign,
    score_mode = score_mode,
    vars_to_regress = vars_to_regress,
    cancer_type = cancer_type,
    pseudobulk = pseudobulk,
    group_by = group_by,
    sample_metadata = sample_metadata,
    percentile = percentile,
    permutation_n = permutation_n,
    permutation_seed = permutation_seed %||% seed,
    verbose = verbose,
    seed = seed,
    ...
  )

  if (verbose) {
    ts_cli$cli_alert_info("Deriving reference from bulk data")
  }

  reference <- PhenoMapR::derive_reference_from_bulk(
    bulk_expression = matched_bulk,
    phenotype = p$pheno_df,
    sample_id_column = "cell_id",
    phenotype_column = p$phenotype_column,
    phenotype_type = p$phenotype_class,
    survival_time = p$survival_time,
    survival_event = p$survival_event,
    normalize = p$normalize,
    platform = p$platform,
    probe_annotation = p$probe_annotation,
    hugo_species = p$hugo_species,
    binary_positive_reference = p$binary_positive_reference,
    verbose = p$verbose
  )

  if (verbose) {
    ts_cli$cli_alert_info("Scoring single cells")
  }

  scores <- PhenoMapR::PhenoMap(
    expression = sc_data,
    reference = reference,
    assay = p$assay,
    layer = p$layer,
    cancer_type = p$cancer_type,
    z_score_cutoff = p$z_score_cutoff,
    pseudobulk = p$pseudobulk,
    group_by = p$group_by,
    # slot = NULL,
    verbose = p$verbose,
    reference_sign = p$reference_sign,
    score_mode = p$score_mode,
    vars_to_regress = p$vars_to_regress,
    sample_metadata = p$sample_metadata,
    cell_id_column = "cell_id",
    permutation_n = p$permutation_n,
    permutation_seed = p$permutation_seed
  )

  if (verbose) {
    ts_cli$cli_alert_info("Defining phenotype groups")
  }

  groups <- PhenoMapR::define_phenotype_groups(
    scores = scores,
    percentile = p$percentile,
    score_columns = NULL
  )

  group_col <- groups[, 2L] |>
    dplyr::recode_values(
      "Most Adverse" ~ "Positive",
      "Most Favorable" ~ "Negative",
      "Other" ~ "Neutral",
      default = NA_character_
    )

  scores$PhenoMapR <- group_col

  sc_data <- SeuratObject::AddMetaData(sc_data, metadata = scores)
  sc_data <- SigBridgeRUtils::AddMisc(
    seurat_obj = sc_data,
    PhenoMapR = S7::props(p$cache_config),
    cover = FALSE
  )

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green("Done"))
  }

  list(
    scRNA_data = sc_data,
    phenomap_scores = scores,
    reference = reference
  )
}
