# ---- 2. Do scPP ----

#' Validate and Prepare DoscPP Parameters
#'
#' @description
#' Internal helper that handles package installation checks, input validation,
#' and default-value resolution for [DoscPP()].
#'
#' @param matched_bulk,sc_data,phenotype,label_type,phenotype_class
#'   Forwarded from [DoscPP()].
#' @param ref_group,Log2FC_cutoff,estimate_cutoff,probs
#'   Forwarded from [DoscPP()].
#' @param ... Additional dots forwarded from [DoscPP()].
#'
#' @return A named list with elements: `phenotype_class`, `verbose`, `parallel`,
#'   `seed`, `assay`.
#'
#' @keywords internal
#' @family scPP
ValidatescPPParams <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type,
  phenotype_class,
  ref_group,
  Log2FC_cutoff,
  estimate_cutoff,
  probs,
  ...
) {
  # -- package checks -------------------------------------------------------
  check_installed("dplyr")
  check_installed("scPAS", action = \(pkg, ...) {
    check_installed("pak")
    pak::pak("Exceret/scPP")
  })

  # -- input validation -----------------------------------------------------
  chk::chk_is(sc_data, "Seurat")
  chk::chk_character(label_type)
  phenotype_class <- tolower(phenotype_class)
  phenotype_class <- arg_match(
    phenotype_class,
    c("binary", "continuous", "survival")
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

  # scPP is more strict than Scissor and scPAS
  if (phenotype_class == "survival") {
    if (!all(rownames(phenotype) == colnames(matched_bulk))) {
      Abort(
        "Please check the rownames of {.var phenotype} and colnames of {.var bulk_dataset},\
         they should be the same."
      )
    }
  } else {
    if (!all(names(phenotype) == colnames(matched_bulk))) {
      Abort(
        "Please check the names of {.var phenotype} and colnames of {.var bulk_dataset},\
         they should be the same."
      )
    }
  }

  # -- process dots ---------------------------------------------------------
  dots <- list2(...)
  verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption("verbose")
  parallel <- (dots$parallel %||% FALSE)
  seed <- dots$seed %||% SigBridgeRUtils::getFuncOption("seed")
  assay <- dots$assay %||% "RNA"
  save_cache <- dots$save_cache
  load_cache <- dots$load_cache

  # -- build cache config ---------------------------------------------------
  cache_config <- ScreenMethodConfig(
    method_name = "scPP",
    param = get_env_vars(exclude = c("matched_bulk", "sc_data", "phenotype")),
    phenotype_class = phenotype_class,
    label_type = label_type
  )

  get_env_vars(exclude = c("matched_bulk", "sc_data", "phenotype"))
}

#' @title Perform scPP screening analysis
#'
#' @description
#' This function performs scPP screening on single-cell data using matched bulk data and phenotype information.
#' It supports binary, continuous, and survival phenotype types.
#'
#'
#' @inheritParams Screen
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
  # -- validate & prepare all parameters -----------------------------------
  p <- exec(ValidatescPPParams, !!!fn_fmls())

  if (p$verbose) {
    ts_cli$cli_alert_info(cli::col_green("Start scPP screening."))
  }

  # decide which type of phenotype data is used
  if (is.vector(phenotype)) {
    Feature <- NULL # suppress checking NOTE
    # The reason why using data.frame instead of vector is to
    # keep the same input and output format with scPP
    phenotype <- as.data.frame(phenotype) |>
      SigBridgeRUtils::Rownames2Col("Sample") |>
      dplyr::rename("Feature" = 2) |>
      dplyr::mutate(Feature = as.numeric(`Feature`))
  }

  if (!is.null(p$load_cache)) {
    cache <- CacheSysCall(
      mode = "load",
      path = p$load_cache,
      cache = p$cache_config,
      verbose = p$verbose,
      timestamp = p$timestamp
    )
    gene_list <- cache$gene_list
    rm(cache)
    if (p$verbose) {
      ts_cli$cli_alert_info("Loaded gene list from cache.")
    }
  } else {
    if (p$verbose) {
      ts_cli$cli_alert_info("Finding overall markers...")
    }

    gene_list <- switch(
      p$phenotype_class,
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
    pos_null <- if ("gene_pos" %chin% names(l)) {
      # Cannot combine the conditions due to the feature of `gene_list`
      if (l[["gene_pos"]] == 0) {
        cli::cli_warn("No significant positive genes found")
        TRUE
      }
      FALSE
    } else {
      FALSE
    }
    neg_null <- if ("gene_neg" %chin% names(l)) {
      if (l[["gene_neg"]] == 0) {
        cli::cli_warn("No significant negative genes found")
        TRUE
      }
      FALSE
    } else {
      FALSE
    }
    if (pos_null && neg_null) {
      cli::cli_warn(
        "scPP is not applicable to the current data. Returning {.val NULL}",
      )
      return(NULL)
    }

    gene_list <- if (is.null(p$save_cache)) {
      # avoid saving cache when loading cache
      cache <- ScreenMethodCache(
        cache_path = p$save_cache,
        cache_config_path = file.path(p$save_cache, "cache_config.json"),
        cache_data = list(phenotype = phenotype, gene_list = gene_list), # NULL is OK; stores the data
        screen_method_config = new_property(class = ScreenMethodConfig)
      )
      CacheSysCall(
        mode = "save",
        cache = cache,
        verbose = p$verbose,
        timestamp = p$dots$timestamp,
      )
    }
  } # end of `is.null(p$load_cache)`

  if (p$verbose) {
    ts_cli$cli_alert_info("Screening")
  }

  # * Start screen
  scPP_result <- ScPP::ScPP.optimized(
    sc_dataset = sc_data,
    geneList = gene_list,
    probs = probs,
    verbose = p$verbose,
    parallel = p$parallel,
    seed = p$seed,
    assay = p$assay
  )
  sc_data[[]] <- scPP_result$metadata
  sc_data <- SigBridgeRUtils::AddMisc(
    sc_data,
    scPP = props(p$cache_config),
    cover = FALSE
  )

  if (p$verbose) {
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
