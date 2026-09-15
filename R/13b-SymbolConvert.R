#' @title Convert Ensembles Version IDs & TCGA Version IDs to Genes in Bulk Expression Data
#'
#' @description
#' Convert Ensembles version IDs and TCGA version IDs to genes. NA values (unmatched IDs)
#' are replaced with the original unconverted IDs (Ensembl version IDs or TCGA version IDs).
#'
#' @param data RNA expression data with rownames as Ensembles version IDs or TCGA version IDs(matrix or data.frame)
#' @param genome_build Genome build of the data. Default: `"hg38"`.
#' @param update_symbol Whether to update gene symbol to latest version. Default: `TRUE`.
#' @param update_symbol_method Method used to update gene symbols, one of `"auto"`,
#'     `"scCustomize"` or `"Seurat"`. `"auto"` prefers `scCustomize` when it is
#'     installed, otherwise falls back to `Seurat`. Only used when
#'     `update_symbol = TRUE`. Default: `"auto"`.
#' @param verbose Whether to print verbose messages. Default: `TRUE`.
#' @param unknown_format `r lifecycle::badge("deprecated")` This argument is deprecated and
#'     ignored. NA values are now replaced with the original unconverted IDs.
#' @param ... No usage
#'
#' @return Data with rownames converted to gene symbols
#'
#' @export
#' @family input_preprocess
#'
SymbolConvert <- function(
  data,
  genome_build = c("hg38", "hg19", "mm10", "mm9", "ce11", "T2T"),
  update_symbol = TRUE,
  update_symbol_method = c("auto", "scCustomize", "Seurat"),
  verbose = SigBridgeRUtils::getFuncOption("verbose"),
  unknown_format = lifecycle::deprecated(),
  ...
) {
  chk::chk_flag(update_symbol)
  chk::chk_flag(verbose)
  dots <- list(...)

  update_symbol_method <- arg_match(
    update_symbol_method,
    c("auto", "scCustomize", "Seurat")
  )

  # -- deprecated `unknown_format` ------------------------------------------
  if (lifecycle::is_present(unknown_format)) {
    lifecycle::deprecate_warn(
      "4.0.0",
      "SymbolConvert(unknown_format = )",
      details = "NA values are now replaced with the original unconverted IDs."
    )
  }

  check_installed("IDConverter")

  genome_build <- arg_match(genome_build)

  is_datalike <- is_2d(data)
  original_ids <- if (is_datalike) {
    rownames(data)
  } else if (is.character(data)) {
    data
  } else {
    Abort(
      "Data must be a matrix-like object with rownames as Ensembles version IDs or TCGA version IDs,\
       or a vector of Ensembles version IDs or TCGA version IDs",
      type = "[DATA ERROR]"
    )
  }
  if (is.null(original_ids)) {
    Abort("Genes are missing in the data", type = "[DATA ERROR]")
  }
  row_names <- if (is_datalike) {
    gsub("\\..*$", "", original_ids)
  } else {
    original_ids
  }

  gene_symbols <- IDConverter::convert_hm_genes(
    row_names,
    genome_build = genome_build
  )

  # * replace NA with original unconverted IDs
  na_count <- if (is_installed("cheapr")) {
    cheapr::na_count(gene_symbols)
  } else {
    sum(is.na(gene_symbols))
  }
  if (na_count > 0L) {
    k <- if (is_installed("cheapr")) {
      cheapr::which_na(gene_symbols)
    } else {
      which(is.na(gene_symbols))
    }
    gene_symbols[k] <- original_ids[k]
    cli::cli_alert_warning(
      "Found {.val {na_count}} NA values (position {.val {k}}) in gene symbols during conversion. Replaced them with the original IDs."
    )
  }

  is_duplicated <- duplicated(gene_symbols)
  if (any(is_duplicated)) {
    where_duplicated <- which(is_duplicated)
    count_duplicated <- sum(is_duplicated)
    gene_symbols[where_duplicated] <- original_ids[where_duplicated]

    cli::cli_alert_warning(
      "Found {.val {count_duplicated}} duplicated gene symbols (position ({.val {where_duplicated}})) during conversion. Replaced them with the original IDs."
    )
  }

  if (update_symbol) {
    # * update to latest symbol
    if (verbose) {
      ts_cli$cli_alert_info("Updating gene symbols to latest version")
    }
    use_scCustomize <- switch(
      update_symbol_method,
      "scCustomize" = {
        check_installed(
          "scCustomize",
          reason = "to update gene symbols with `update_symbol_method = \"scCustomize\"`."
        )
        TRUE
      },
      "Seurat" = FALSE,
      is_installed("scCustomize")
    )
    gene_symbols <- if (use_scCustomize) {
      if (grepl("hg|HG", genome_build)) {
        scCustomize::Updated_HGNC_Symbols(
          gene_symbols,
          case_check_as_warn = TRUE,
          verbose = verbose
        )$Output_Features
      } else {
        scCustomize::Updated_MGI_Symbols(
          gene_symbols,
          verbose = verbose
        )$Output_Features
      }
    } else {
      Seurat::UpdateSymbolList(gene_symbols, verbose = verbose)
    }
  }

  if (is_datalike) {
    rownames(data) <- gene_symbols
    return(data)
  }
  gene_symbols
}
