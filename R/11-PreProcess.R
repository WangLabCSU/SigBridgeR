# ---- Input Data Preprocess Functions ----

#' @title Preprocess single-cell RNA data
#'
#' @description
#' Preprocess single-cell RNA data: Create a Seurat object, normalize data, find variable features, scale data, run PCA, find clusters and neighbors, run tSNE and UMAP.
#'
#' @param sc raw scRNA data read from `anndata::read.h5ad`, contains `$X` matrix and `$obs` metadata
#' @param project The project name of the Seurat object, transfered to `SeuratObject@project`.
#' @param min_cells The minimum number of cells that a feature must be observed in to be included in `VariableFeatures`.
#' @param min_features The minimum number of cells that a feature must be expressed in to be included in `VariableFeatures`.
#' @param normalization_method The method used to normalize the data. Options are `"LogNormalize"` (default), `"CLR"`, `"RC"` or `"RCLogNormalize"`.
#' @param scale_factor A numeric value by which the data will be multiplied. Default is `10000`.
#' @param selection_method The method used to select variable features, transferred to `SeuratObject@assays$RNA@var.features`. Options are `"vst"` (default), `"mean.var.plot"`, `"dispersion"`, `"mean.var.plot"` or `"dispersion.scaled"`.
#' @param resolution The resolution parameter for clustering. Default is `0.6`.
#' @param dims The dimensions to use for clustering. Default is `1:10`.
#' @param verbose Whether to show detailed information. Default is `TRUE`.
#' @param future_global_maxsize The memory allocated.
#' @return A Seurat object
#'
#' @export
SCPreProcess = function(
  sc,
  project = "Scissor_Single_Cell",
  min_cells = 400,
  min_features = 0,
  normalization_method = "LogNormalize",
  scale_factor = 10000,
  selection_method = "vst",
  resolution = 0.6,
  dims = 1:10,
  verbose = TRUE,
  future_global_maxsize = 6 * 1024^3
) {
  library(dplyr)
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package Seurat not installed")
  }
  if (is.null(sc$X)) {
    stop("Input must contain $X matrix")
  }
  if (is.null(cnv_status) || length(cnv_status) > 1) {
    cnv_status = cnv_status[[1]]
  }
  sc_matrix <- if (inherits(sc$X, "sparseMatrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Matrix package required for sparse matrix support")
    }
    Matrix::t(sc$X)
  } else {
    t(sc$X)
  }
  options(future.globals.maxSize = future_global_maxsize) # default 6GB

  sc_seurat = SeuratObject::CreateSeuratObject(
    counts = sc_matrix,
    project = project,
    min.cells = min_cells,
    min.features = min_features
  ) %>%
    Seurat::NormalizeData(
      normalization.method = normalization_method,
      scale.factor = scale_factor,
      verbose = verbose
    ) %>%
    Seurat::FindVariableFeatures(
      selection.method = selection_method,
      verbose = verbose
    ) %>%
    Seurat::ScaleData(verbose = verbose) %>%
    Seurat::RunPCA(
      features = SeuratObject::VariableFeatures(.),
      verbose = verbose
    ) %>%
    Seurat::AddMetaData(metadata = sc$obs)

  n_pcs <- ncol(sc_seurat@reductions$pca)
  if (is.null(dims) || max(dims) > n_pcs) {
    dims <- 1:min(max(dims), n_pcs)
    cli::cli_alert_warning(crayon::yellow(
      "The input dimension is greater than the dimension in PCA. It is now set to the maximum dimension in PCA."
    ))
  }

  process_sc_object <- function(obj) {
    obj <- obj %>%
      Seurat::FindNeighbors(dims = dims, verbose = verbose) %>%
      Seurat::FindClusters(resolution = resolution, verbose = verbose) %>%
      Seurat::RunTSNE(dims = dims) %>%
      Seurat::RunUMAP(dims = dims, verbose = verbose)
    return(obj)
  }

  result = process_sc_object(sc_seurat)

  return(result)
}


#' @title Preprocess bulk expression data
#'
#' @description
#' Preprocess bulk expression data: convert Ensembles version IDs and TCGA version IDs to genes.
#' @param data raw bulk expression data
#'
#' @export
#'
BulkPreProcess = function(data) {
  rownames(data) <- substr(rownames(data), 1, 15)
  options(
    IDConverter.datapath = system.file("extdata", package = "IDConverter")
  )
  rownames(data) <- IDConverter::convert_hm_genes(rownames(data))
  return(data)
}


#' @title Preprocess Mutational Signature Data
#'
#' @description
#' Filters and cleans mutational signature data based on tumor type, column validity thresholds,
#' and signature accuracy. Designed for genomic datasets with columns representing
#' mutational signatures (e.g., SBS, DBS, CNV).
#'
#' @param ms_signature A data frame containing raw mutational signature data.
#'    Must include:
#'    - A column matching tumor type (auto-detected)
#'    - A column with signature accuracy scores (auto-detected)
#'    - Numeric columns for signatures (filtered by `col_thresh` and `ms_search_pattern`)
#' @param filter_tumor_type Character. Specify tumor type(s) to retain.
#'    Special values:
#'    - `"all"` (default): Keeps all tumor types
#'    - `"luad"`: Keeps tumor types containing the characters "luad" (case-insensitive)
#' @param col_thresh Numeric [0-1]. Minimum proportion of non-zero values required to
#'    retain a signature column (default: 0.05). Columns with >95% zeros will be dropped.
#' @param accuracy_thresh Numeric [0-1]. Minimum accuracy score to retain samples
#'    (default: 0, keeping all samples regardless of accuracy).
#' @param ms_search_pattern Regex pattern to identify mutational signature columns
#'    (default: "SBS|DBS|CN|CNV|SV|ID|INDEL", case-insensitive).
#'
#' @return A filtered data frame with:
#'    - Rows: Samples passing tumor type and accuracy filters
#'    - Columns: Annotations + signature columns passing `col_thresh` and pattern match
#'    Returns `NULL` with warning if no signature columns remain.
#'
#' @section Processing Steps:
#' 1. **Column Filtering**:
#'    - Keeps annotation columns (character/non-numeric)
#'    - Drops numeric columns with >(1-`col_thresh`) zeros
#' 2. **Row Filtering**:
#'    - Subsets by `filter_tumor_type`
#'    - Drops samples with accuracy < `accuracy_thresh`
#' 3. **Validation**:
#'    - Checks if any signature columns remain
#'    - Warns if all filtered out
#'
#' @export
#' @importFrom dplyr filter select
#' @importFrom cli cli_alert_info cli_alert_warning
#' @importFrom crayon green red bold yellow black
#' @importFrom glue glue
#'
MSPreProcess <- function(
  ms_signature,
  filter_tumor_type = "all",
  col_thresh = 0.05,
  accuracy_thresh = 0,
  ms_search_pattern = "SBS|DBS|CN|CNV|SV|ID|INDEL"
) {
  library(dplyr)

  keep_colnames <- lapply(X = names(ms_signature), FUN = function(col_name) {
    if (grepl(tolower(ms_search_pattern), tolower(col_name))) {
      cli::cli_alert_info(
        "{crayon::red(col_name)} not kept, because it does not match the pattern"
      )
      return(NULL)
    }
    col <- ms_signature[[col_name]]

    if (is.character(col)) {
      cli::cli_alert_info("{crayon::green(col_name)} kept")
      return(col_name) # keep annotational columns
    } else if (is.numeric(col)) {
      ifelse(
        # Retain cols where the proportion of non-zero data in the ms exceeds the threshold.
        mean(col != 0) >= col_thresh, # equal to: sum(col != 0)/length(col) >= col_thresh
        {
          cli::cli_alert_info("{crayon::green(col_name)} kept")
          return(col_name)
        },
        {
          cli::cli_alert_info(
            "{crayon::red(col_name)} not kept, because the valid data in the column < {crayon::bold(col_thresh)}"
          )
          return(NULL)
        }
      )
    } else {
      return(NULL)
    }
  })
  # *filter accracy
  accuracy_column <- grep("[aA][cC]{2}.*", colnames(ms_signature), value = TRUE)
  # *filter cancer type
  tumor_type_col <- grep(
    "[tT]umor|[Cc]ancer",
    colnames(ms_signature),
    value = TRUE
  )

  if (tolower(filter_tumor_type) %in% c("all", "all tumor", "all types")) {
    filter_tumor_type = unique(ms_signature[[tumor_type_col]])
  }

  processed_ms_signature <- ms_signature %>%
    dplyr::filter(
      .data[[accuracy_column]] >= accuracy_thresh,
      grepl(
        pattern = glue::glue(tolower(filter_tumor_type), .sep = "|"),
        x = tolower(.data[[tumor_type_col]])
      )
    ) %>%
    dplyr::select(unlist(keep_colnames))

  # final check for data availability
  if (sum(grepl(ms_search_pattern, colnames(processed_ms_signature))) == 0) {
    cli::cli_alert_warning(crayon::yellow(glue::glue(
      "{crayon::bold('All data columns were filtered out, possible reasons:')}",
      "1. Threshold too high (current: {crayon::black(col_thresh)})",
      "2. No matching data columns existed (search pattern: {crayon::black(ms_search_pattern)})",
      "3. All values were zero",
      .sep = "\n"
    )))
  }

  return(processed_ms_signature)
}


#' @description
#' There are differences in TCGA samples across different datasets. This function is used to find common samples for inner join, and other samples will not be retained.
#' @param col_id description
MatchSample = function(
  ms_signature,
  TCGA_exp_count,
  col_id,
  ms_status_thresh = 0
) {
  require(dplyr)
  stopifnot(col_id %in% seq_along(ms_signature))

  # find sample column
  sample_colname = grep(
    "sample",
    colnames(ms_signature),
    ignore.case = TRUE,
    value = TRUE
  )
  # minimum length can be used to match samples
  trunc_length <- min(
    min(nchar(ms_signature[[sample_colname]])),
    min(nchar(colnames(TCGA_exp_count))),
    15
  )
  # convert interested column to binary variable
  processed_ms_signature <- ms_signature %>%
    dplyr::mutate(
      # Samples exceeding the threshold will be retained and set logical
      !!sample_colname := substr(.[[sample_colname]], 1, trunc_length),
      ms_status = as.integer(.[[col_id]] > ms_status_thresh)
    )
  colnames(TCGA_exp_count) = substr(colnames(TCGA_exp_count), 1, trunc_length)
  cli::cli_alert_info("Sample match length: {trunc_length}")
  # find common sample names
  cm_samples = intersect(
    processed_ms_signature[[sample_colname]],
    colnames(TCGA_exp_count)
  )
  cli::cli_alert_info("Sample match: {length(cm_samples)} common samples")

  if (length(cm_samples) == 0) {
    stop("No common sample found in data")
  }

  match_result = list(
    phenotype = processed_ms_signature %>%
      dplyr::filter(.[[sample_colname]] %in% cm_samples) %>%
      dplyr::arrange(factor(.[[sample_colname]], levels = cm_samples)) %>%
      {
        stats::setNames(.$ms_status, .[[sample_colname]])
      },
    matched_TCGA_exp_count = TCGA_exp_count[, cm_samples, drop = FALSE],
    ms_select = names(ms_signature)[[col_id]]
  )

  return(match_result)
}

# ------------------- other function ----------------------

IdentifyDataColumn <- function(data, thresh = 0.5) {
  pattern_based <- sapply(data, function(x) {
    if (is.character(x)) {
      numeric_pat <- grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", x)
      mean(numeric_pat, na.rm = TRUE) > 0.5
    } else {
      TRUE
    }
  })
  return(pattern_based)
}
