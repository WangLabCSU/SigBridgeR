#' @title Perform scPP screening analysis
#'
#' @description
#' This function performs scPP screening on single-cell data using matched bulk data and phenotype information.
#' It supports binary, continuous, and survival phenotype types.
#'
#' @usage
#' DoscPP(
#'   matched_bulk,
#'   sc_data,
#'   phenotype,
#'   label_type = "scPP",
#'   phenotype_class = c("binary", "continuous", "survival"),
#'   ref_group = 0,
#'   Log2FC_cutoff = 0.585,
#'   estimate_cutoff = 0.2,
#'   probs = c(0.2, NULL),
#'   ...
#' )
#'
#' @param matched_bulk Bulk expression data (genes × samples) where:
#'        - Column names must match `phenotype` row names
#' @param sc_data Seurat object containing preprocessed single-cell data:
#'        - Normalized counts in `RNA` assay
#' @param phenotype Data frame or tibble or named vector with:
#'        - Rownames matching `matched_bulk` columns
#'        - For survival: must contain time and status columns
#' @param label_type Character specifying phenotype label type (e.g., "SBS1"), stored in `scRNA_data@misc`
#' @param phenotype_class Analysis type (case-sensitive):
#'        - `"Binary"`: Case-control studies (e.g., tumor/normal)
#'        - `"Continuous"`: Quantitative traits (e.g., drug response)
#'        - `"Survival"`: Time-to-event data (requires time/status columns)
#' @param ref_group Reference group or baseline for **binary** comparisons, e.g. "Normal" for Tumor/Normal studies and 0 for 0/1 case-control studies. (default: 0)
#' @param Log2FC_cutoff Minimum log2 fold-change for binary markers (default: 0.585)
#' @param estimate_cutoff Effect size threshold for **continuous** traits (default: 0.2)
#' @param probs A numeric value indicating the quantile cutoff for cell classification. This parameter can also be a numeric vector, in which case an optimal threshold will be selected based on the AUC and enrichment score.(default: `0.2`)
#' @param ... Additional arguments. Currently supports:
#'    - `verbose`: Logical indicating whether to print progress messages. Defaults to `TRUE`.
#'    - `seed`: For reproducibility, default is `123L`
#'    - `assay`: The assay to use for single-cell data. Defaults to `RNA`.
#'
#'
#' @return A list containing:
#' \describe{
#'   \item{scRNA_data}{Seurat object with added metadata:
#'     \describe{
#'       \item{ScPP}{"Positive"/"Negative"/"Neutral" classification}
#'     }
#'   }
#'   \item{gene_list}{List of genes used for screening}
#' }
#'
#' @section Algorithm Steps:
#' 1. Data Validation: Checks sample alignment between bulk and phenotype data
#' 2. Marker Selection: Identifies phenotype-associated genes from bulk data
#' 3. Single-cell Screening: Projects bulk markers onto single-cell data
#' 4. Cell Classification: Categorizes cells based on phenotype association
#'
#' @section Reference:
#' WangX-Lab/ScPP \[Internet\]. \[cited 2025 Aug 31\]. Available from: https://github.com/WangX-Lab/ScPP
#'
#' @examples
#' \dontrun{
#' # Binary phenotype analysis
#' res <- DoscPP(
#'   matched_bulk = bulk_data,
#'   sc_data = seurat_obj,
#'   phenotype = ms_data,
#'   label_type = "SBS1",
#'   phenotype_class = "Binary"
#' )
#'
#' # Survival analysis
#' surv_res <- DoscPP(
#'   sc_data = seurat_obj,
#'   matched_bulk = bulk_data,
#'   phenotype = surv_df,
#'   label_type = "OS_status",
#'   phenotype_class = "Survival"
#' )
#' }
#'
#' @export
#' @family screen_method
#' @family scPP
#'
DoscPP <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = "scPP",
  phenotype_class = c("binary", "continuous", "survival"),
  ref_group = 0L,
  Log2FC_cutoff = 0.585,
  estimate_cutoff = 0.2,
  probs = c(0.2, NULL),
  ...
) {
  CheckInstalled("Exceret/ScPP")
  chk::chk_is(sc_data, "Seurat")
  chk::chk_character(label_type)
  phenotype_class <- tolower(phenotype_class)
  phenotype_class <- SigBridgeRUtils::MatchArg(
    phenotype_class,
    c("binary", "continuous", "survival"),
    NULL
  )
  chk::chk_number(ref_group)
  chk::chk_range(Log2FC_cutoff)
  chk::chk_range(estimate_cutoff)
  if (!is.null(probs)) {
    chk::chk_range(probs, range = c(0, 0.5))
  }
  # scPP can't tolerate NA
  chk::chk_not_any_na(matched_bulk)
  chk::chk_not_any_na(phenotype)

  # robust, scPP is more strict than scissor and scPAS
  if (phenotype_class == "survival") {
    if (!all(rownames(phenotype) == colnames(matched_bulk))) {
      cli::cli_abort(c(
        "x" = "Please check the rownames of {.var phenotype} and colnames of {.var bulk_dataset}, they should be the same."
      ))
    }
  } else {
    if (!all(names(phenotype) == colnames(matched_bulk))) {
      cli::cli_abort(c(
        "x" = "Please check the names of {.var phenotype} and colnames of {.var bulk_dataset}, they should be the same."
      ))
    }
  }

  # * default options
  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption("verbose")
  parallel <- (dots$parallel %||% FALSE) &
    !inherits(future::plan("list")[[1]], "sequential")
  seed <- dots$seed %||% SigBridgeRUtils::getFuncOption("seed")
  assay <- dots$assay %||% "RNA"

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green("Start scPP screening."))
    ts_cli$cli_alert_info("Finding overall markers...")
  }

  # decide which type of phenotype data is used
  if (is.vector(phenotype)) {
    # The reason why using data.frame instead of vector is to
    # keep the same input and output format with scPP
    phenotype <- as.data.frame(phenotype) %>%
      SigBridgeRUtils::Rownames2Col("Sample") %>%
      dplyr::rename("Feature" = 2) %>%
      dplyr::mutate(Feature = as.numeric(`Feature`))
  }

  gene_list <- switch(
    phenotype_class,
    "binary" = {
      if (estimate_cutoff != 0.2) {
        cli::cli_warn(
          "The parameters {.arg estimate_cutoff} are not used for survival analysis. Ignore it"
        )
      }
      ScPP::marker_Binary.optimized(
        bulk_data = matched_bulk,
        features = phenotype,
        ref_group = ref_group,
        Log2FC_cutoff = Log2FC_cutoff
      )
    },
    "continuous" = {
      if (Log2FC_cutoff != 0.585) {
        cli::cli_warn(
          "The parameters {.arg Log2FC_cutoff} are not used for survival analysis. Ignore it"
        )
      }
      ScPP::marker_Continuous.optimized(
        bulk_data = matched_bulk,
        features = phenotype$Feature,
        method = "spearman",
        estimate_cutoff = estimate_cutoff
      )
    },
    "survival" = {
      if (estimate_cutoff != 0.2) {
        cli::cli_warn(
          "The parameters {.arg estimate_cutoff} are not used for survival analysis. Ignore it"
        )
      }
      if (Log2FC_cutoff != 0.585) {
        cli::cli_warn(
          "The parameters {.arg Log2FC_cutoff} are not used for survival analysis. Ignore it"
        )
      }
      ScPP::marker_Survival2(
        bulk_data = matched_bulk,
        survival_data = phenotype
      )
    }
  )

  l <- lapply(gene_list, length)
  pos_null <- FALSE
  neg_null <- FALSE
  if ("gene_pos" %chin% names(l)) {
    # Cannot combine the conditions due to the feature of `gene_list`
    if (l[["gene_pos"]] == 0 && verbose) {
      ts_cli$cli_alert_info("No significant positive genes found")
      pos_null <- TRUE
    }
  }
  if ("gene_neg" %chin% names(l)) {
    if (l[["gene_neg"]] == 0 && verbose) {
      ts_cli$cli_alert_info("No significant negative genes found")
      neg_null <- TRUE
    }
  }
  if (pos_null && neg_null) {
    cli::cli_warn(
      "scPP is not applicable to the current data. Returning {.val NULL}",
    )
    return(NULL)
  }

  if (verbose) {
    ts_cli$cli_alert_info("Screening")
  }

  # * Start screen
  scPP_result <- ScPP::ScPP.optimized(
    sc_dataset = sc_data,
    geneList = gene_list,
    probs = probs,
    verbose = verbose,
    parallel = parallel,
    seed = seed,
    assay = assay
  )
  sc_data[[]] <- scPP_result$metadata
  sc_data <- SigBridgeRUtils::AddMisc(
    sc_data,
    scPP_type = label_type,
    scPP_para = list(
      ref_group = ref_group,
      Log2FC_cutoff = Log2FC_cutoff,
      estimate_cutoff = estimate_cutoff,
      probs = probs
    ),
    cover = FALSE
  )

  if (verbose) {
    ts_cli$cli_alert_success(cli::col_green("scPP screening done."))
  }

  list(
    scRNA_data = sc_data,
    gene_list = list(
      genes_pos = scPP_result$Genes_pos,
      genes_neg = scPP_result$Genes_neg
    )
  )
}
