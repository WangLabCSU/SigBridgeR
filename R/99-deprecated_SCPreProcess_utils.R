# nocov start

#' @keywords internal
#' @family single_cell_preprocess
has_pattern <- function(qc_list) {
  if (!is.list(qc_list)) {
    return(FALSE)
  }
  !is.null(qc_list) && !is.null(qc_list$pattern) && length(qc_list$pattern) > 0
}

#' @keywords internal
#' @family single_cell_preprocess
is_filtering <- function(filter_list) {
  if (!is.list(filter_list)) {
    return(FALSE)
  }
  !is.null(filter_list)
}

#' @keywords internal
#' @family single_cell_preprocess
modifyList2 <- function(x, val) {
  chk::chk_list(x)
  chk::chk_list(val)

  purrr::imap(val, function(v_val, key) {
    if (key %chin% names(x) && is.list(x[[key]]) && is.list(v_val)) {
      return(utils::modifyList(x[[key]], v_val))
    }
    v_val
  })
}

#' @keywords internal
#' @family single_cell_preprocess
handle_usr_params <- function(usr_params) {
  default <- list(
    # * CreateSeuratObject
    o = list(
      project = "SC_Screen_Proj",
      min.cells = 400L
    ),
    # * NormalizeData
    n = list(),
    # * ScaleData
    s = list(),
    # * FindVariableFeatures
    v = list(),
    # * RunPCA
    p = list(),
    # * FindNeightbors
    e = list(),
    # * FindClusters
    c = list(
      resolution = 0.6
    ),
    # * RunTSNE
    t = list(),
    # * RunUMAP
    u = list()
    # * SCTransform
    # r = list()
  )
  modifyList2(default, usr_params)
}


#' @keywords internal
#' @family single_cell_preprocess
GetVars2Regress <- function(seurat_obj, verbose = TRUE) {
  if (!"qc_colnames" %chin% names(seurat_obj@misc)) {
    cli::cli_warn(
      "[{.fun GetVars2Regress}]: No QC columns found for regression in {.fun SCTransform}. Please run {.fun QCPatternDetect} first. Returning {.val NULL}."
    )

    return(NULL)
  }

  qc_cols <- unlist(seurat_obj@misc$qc_colnames, use.names = FALSE)
  meta <- seurat_obj[[]]

  existing_cols <- qc_cols[qc_cols %chin% colnames(meta)]

  if (length(existing_cols) == 0) {
    cli::cli_warn(
      "[{.fun GetVars2Regress}]: No matching QC columns found in seurat@meta.data Returning {.val NULL}."
    )

    return(NULL)
  }

  # remove constant columns (variance = 0)
  valid_cols <- existing_cols[vapply(
    existing_cols,
    function(col) {
      x <- meta[[col]]
      !all(is.na(x)) && stats::var(x, na.rm = TRUE) > 0
    },
    logical(1L)
  )]

  if (length(valid_cols) == 0) {
    cli::cli_warn(
      "[{.fun GetVars2Regress}]: All QC columns are constant, skipping regression. Returning {.val NULL}."
    )
    return(NULL)
  }

  if (verbose && length(valid_cols) > 0) {
    cli::cli_text(
      "Using QC columns {valid_cols} for regression in {.fun SCTransform}"
    )
  }

  valid_cols
}

#' @title Filter tumor cells (internal)
#'
#' @description
#' An internal function that filters tumor cells from a Seurat object based on
#' metadata column values. This function identifies tumor cells using pattern
#' matching on cell type labels and creates a subset containing only tumor cells.
#' It also records dimension information before and after filtering for
#' traceability.
#'
#' @param obj Seurat object with a column to filter out tumor cells.
#' @param column2only_tumor Name of the column to filter out tumor cells.
#' @param verbose logical, whether to print progress messages
#'
#' @return A Seurat object containing only tumor cells, with the following
#'         attributes stored in `@misc`:
#'
#'   * `self_dim`: Dimensions of the filtered object
#'   * `raw_dim`: Original dimensions before filtering
#'   * `column2only_tumor`: The column name used for filtering
#'
#'   If `column2only_tumor` is `NULL` or the specified column is not found,
#'   returns the original object unchanged.
#'
#' @keywords internal
#' @family single_cell_preprocess
#'
FilterTumorCell <- function(
  obj,
  column2only_tumor = NULL,
  verbose = TRUE
) {
  raw_dim <- dim(obj)
  obj <- SigBridgeRUtils::AddMisc(obj, self_dim = raw_dim, cover = TRUE)

  if (is.null(column2only_tumor)) {
    return(obj)
  }
  if (!column2only_tumor %chin% colnames(obj[[]])) {
    cli::cli_warn(
      "Column '{.emph {column2only_tumor}}' not found, skip tumor cell filtering"
    )
    return(obj)
  }
  if (verbose) {
    cli::cli_alert_info(
      "Filtering tumor cells with '{.emph {column2only_tumor}}'..."
    )
  }

  labels <- obj[[column2only_tumor]][[1]]
  tumor_cells <- grepl(
    "^[Tt]umo.?r|[Cc]ancer[Mm]alignant|[Nn]eoplasm|[Tt]um|1",
    labels
  )

  obj <- obj[, tumor_cells]

  SigBridgeRUtils::AddMisc(
    seurat_obj = obj,
    raw_dim = raw_dim,
    self_dim = dim(obj),
    column2only_tumor = column2only_tumor,
    cover = TRUE
  )
}

#' @title Calculate Percentage of Features Matching Patterns
#'
#' @description
#' This function calculates the percentage of counts coming from features matching
#' specified patterns (e.g., mitochondrial genes, ribosomal genes) and adds them
#' as metadata columns to the Seurat object.
#'
#' @param obj A seurat object.
#' @param pattern A character vector or list containing regex patterns to identify mitochondrial
#'    genes, ribosomal protein genes, or other unwanted genes, as well as combinations
#'    of these genes. Customized patterns are supported.
#' @param verbose logical, whether to print progress messages
#' @param ... Additional arguments passed to \code{\link[Seurat]{PercentageFeatureSet}}
#'
#' @details
#' The function automatically generates friendly column names based on the patterns:
#' - "mt" for mitochondrial patterns
#' - "rp" for ribosomal patterns
#' - "rrna" for ribosomal RNA patterns
#' - For combined patterns (using |), creates names like "mt_rp"
#' - For other patterns, creates cleaned lowercase names
#'
#' @family single_cell_preprocess
#' @export
#'
QCPatternDetect <- function(
  obj,
  pattern = c("^MT-", "^mt-", "^RP[SL]", "^MT-|^RP[SL]"),
  verbose = TRUE,
  ...
) {
  if (verbose) {
    cli::cli_text(
      "Using QC patterns {.arg {pattern}} to detect metrics"
    )
  }

  # if `pattern` is a list, convert it to a character vector
  patterns <- unlist(pattern)
  colname_mapping <- stats::setNames(
    paste0("percent.", purrr::map_chr(patterns, Pattern2Colname)),
    patterns
  )
  # Maybe these filter already exist in the object
  existing_cols <- colnames(obj[[]])
  existing <- colname_mapping %chin% existing_cols
  patterns_to_process <- names(colname_mapping)[!existing]

  if (verbose && any(existing)) {
    skipped_cols <- colname_mapping[existing]
    skipped_patterns <- names(skipped_cols)

    purrr::walk2(
      skipped_cols,
      skipped_patterns,
      ~ cli::cli_warn(
        "[{.fun QCPatternDetect}]: Column {.val {.x}} already exists. Skipping pattern: {.val {.y}}"
      )
    )
  }

  obj <- purrr::reduce(
    .x = patterns_to_process,
    .f = function(obj_acc, pat) {
      col_name <- colname_mapping[[pat]]
      obj_acc[[col_name]] <- Seurat::PercentageFeatureSet(
        obj_acc,
        pattern = pat,
        ...
      )
      obj_acc
    },
    .init = obj
  )

  # Record these colnames to misc slot for further data filter
  obj@misc$qc_colnames <- unname(colname_mapping)

  obj
}

#' @title Filter Seurat object cells by QC metrics
#'
#' @description
#' Filters cells based on nFeature_RNA and optional QC metrics (e.g. percent.mt, percent.rp),
#' defined in `seurat_obj@misc$qc_colnames` (See \code{\link{QCPatternDetect}}). Only metrics with non-constant, non-all-zero values are used.
#'
#' @param seurat_obj A \code{Seurat} object.
#' @param data_filter.thresh A named list with thresholds. Default:
#'   \code{list(
#'     nFeature_thresh = c(200L, 6000L),
#'     nCount_thresh    = c(500L, 25000L),
#'     percent.mt = 20L,
#'     percent.rp = 60L
#'   )}.
#' @param verbose Logical; whether to print progress messages via \code{cli}.
#' @param ... No use
#'
#' @return A filtered \code{Seurat} object.
#' @export
#'
#'
QCFilter <- function(
  seurat_obj,
  data_filter.thresh = list(
    assay = SeuratObject::DefaultAssay(seurat_obj),
    nFeature_thresh = c(200L, 6000L),
    nCount_thresh = c(500L, 50000L),
    percent.mt = 20L,
    percent.rp = 60L
  ),
  verbose = TRUE,
  ...
) {
  check_installed("dplyr")
  chk::chk_is(seurat_obj, "Seurat")
  chk::chk_list(data_filter.thresh)

  defaults <- list(
    assay = names(seurat_obj@assays)[[1]],
    nFeature_thresh = c(200L, 6000L),
    nCount_thresh = c(500L, 50000L),
    percent.mt = 20L,
    percent.rp = 60L
  )
  thresh <- utils::modifyList(defaults, data_filter.thresh)

  if (verbose) {
    feature_name <- paste0("nFeature_", thresh$assay)
    count_name <- paste0("nCount_", thresh$assay)
    cli::cli_text(
      "Filtering cells by {.arg {feature_name}}, {.arg {count_name}} and {.arg QC metrics}"
    )
  }

  # * Ensure nFeature_thresh is length-2 integer
  chk::chk_length(thresh$nFeature_thresh, 2)
  chk::chk_length(thresh$nCount_thresh, 2)

  thresh$nFeature_thresh <- as.integer(thresh$nFeature_thresh)
  thresh$nCount_thresh <- as.integer(thresh$nCount_thresh)

  # filter expr is a string of the form
  # which is used to subset the Seurat object
  nfeat_condition <- expr(
    !!as.symbol(paste0("nFeature_", thresh$assay)) >
      !!thresh$nFeature_thresh[1] &
      !!as.symbol(paste0("nFeature_", thresh$assay)) <
        !!thresh$nFeature_thresh[2]
  )
  ncount_condition <- expr(
    !!as.symbol(paste0("nCount_", thresh$assay)) > !!thresh$nCount_thresh[1] &
      !!as.symbol(paste0("nCount_", thresh$assay)) < !!thresh$nCount_thresh[2]
  )

  # see `QCPatternDetect()` for the column names generation
  qc_colnames <- NULL
  if ("qc_colnames" %chin% names(seurat_obj@misc)) {
    qc_colnames <- unlist(seurat_obj@misc$qc_colnames, use.names = FALSE)
  }
  if (is.null(qc_colnames) || length(qc_colnames) == 0) {
    qc_colnames <- character(0)
  }

  # --- Validate QC columns and build conditions -------------------------------
  get_qc_condition <- function(qc_colnames, meta, thresh) {
    if (length(qc_colnames) == 0) {
      return(NULL)
    }

    present_cols <- qc_colnames[qc_colnames %chin% colnames(meta)]

    valid_cols <- present_cols[vapply(
      X = present_cols,
      FUN = function(col) {
        x <- meta[[col]]
        has_vals <- any(!is.na(x) & x >= 0, na.rm = TRUE)
        non_0_var <- if (is_installed("cheapr")) {
          cheapr::var_(x, na.rm = TRUE) > 0
        } else {
          stats::var(x, na.rm = TRUE) > 0
        }

        has_vals & non_0_var
      },
      FUN.VALUE = logical(1L)
    )]

    if (length(valid_cols) == 0) {
      return(NULL)
    }

    purrr::map(valid_cols, function(col) {
      if (!col %chin% names(thresh)) {
        cli::cli_warn(
          "[{.fun QCFilter}]: No threshold provided for {.val {col}}; skipping."
        )
        return(NULL)
      }
      expr(!!dplyr::sym(col) < !!thresh[[col]])
    }) |>
      purrr::compact()
  }
  # be aware of expr is not supported in `subset()`
  meta <- seurat_obj[[]]
  qc_conds <- get_qc_condition(qc_colnames, meta, thresh)

  all_conds <- c(list(nfeat_condition, ncount_condition), qc_conds)
  if (length(all_conds) == 0) {
    Abort(
      "[{.fun QCFilter}]: No valid filtering conditions generated."
    )
  }

  full_expr <- purrr::reduce(
    all_conds,
    .f = function(x, y) expr(!!x & !!y)
  )

  # Evaluate safely
  logical_vec <- base::with(data = meta, expr = eval_tidy(full_expr))

  if (!is.logical(logical_vec) || length(logical_vec) != nrow(meta)) {
    Abort(
      "[{.fun QCFilter}]: filtering condition did not produce a logical vector of length {.val {nrow(meta)}}."
    )
  }

  keep_cells <- rownames(meta)[logical_vec]

  if (verbose) {
    n_kept <- length(keep_cells)
    n_total <- nrow(meta)
    pct_off <- if (n_total > 0) 100 * (1 - n_kept / n_total) else 0
    cli::cli_text(sprintf(
      "Kept  %d/%d (%.2f%% off) cells after filtering",
      n_kept,
      n_total,
      pct_off
    ))
  }

  subset(seurat_obj, cells = keep_cells)
}

#' @title convert regex patterns to column names (internal)
#'
#' @description
#' An internal utility function that converts regular expression patterns used
#' for quality control in single-cell RNA-seq analysis into standardized column
#' names for Seurat object metadata. This function handles both single patterns
#' and combined patterns with logical OR operators.
#'
#' @param pat A character string containing a regular expression pattern or
#'            multiple patterns combined with `|` (OR operator).
#'
#' @return A character string with the standardized column name derived from
#'         the input pattern(s). The output follows these conventions:
#'         - Lowercase letters only
#'         - Special characters removed or replaced with underscores
#'         - Common patterns mapped to standardized abbreviations
#'         - Combined patterns sorted alphabetically and joined with underscores
#'
#' @examples
#' \dontrun{
#' # Internal usage examples
#' Pattern2Colname("^MT-")                    # Returns "mt"
#' Pattern2Colname("^RP[LS]")                 # Returns "rp"
#' Pattern2Colname("^[rt]rna")                # Returns "rt_rna"
#' Pattern2Colname("^MT-|^RP[LS]")            # Returns "mt_rp"
#' Pattern2Colname("^HB[AB]?")                # Returns "hb_ab"
#' Pattern2Colname("Custom_Pattern[0-9]+")    # Returns "custom_pattern_0_9"
#' }
#'
#' @export
#' @family single_cell_preprocess
#'
Pattern2Colname <- function(pat) {
  pat_lower <- tolower(pat)

  if (grepl("\\|", pat_lower)) {
    # Handle combined patterns (with | separator)
    parts <- strsplit(pat_lower, "|", fixed = TRUE)[[1]]
    names <- purrr::map_chr(parts, function(p) {
      dplyr::case_when(
        grepl("mt", p) ~ "mt",
        grepl("rp", p) ~ "rp",
        TRUE ~ tolower(gsub("[^[:alnum:]]", "", p))
      )
    })

    return(paste(sort(unique(names)), collapse = "_"))
  }
  # Handle single patterns
  dplyr::case_when(
    grepl("mt", pat_lower) ~ "mt",
    grepl("rp", pat_lower) ~ "rp",
    TRUE ~ {
      clean <- gsub("[^[:alnum:]]", "_", pat_lower)
      clean <- gsub("_+", "_", clean) # Collapse multiple underscores
      clean <- gsub("^_+|_+$", "", clean) # Trim leading/trailing underscores
      tolower(clean)
    }
  )
}
