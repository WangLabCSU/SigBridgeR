# ---- Input Data Preprocess Functions ----

#' @description
#' Preprocess single-cell RNA data: Create a Seurat object, filter out the tumor cells from it, and then create another Seurat object using these tumor cells.
#' @param sc raw scRNA data
#' @param cnv_status The criterion for filtering tumor cells from the Seurat object, i.e. `cnv_status == "tumor"`. Additionally, another option "normal" is provided (although no use).
#' @param future_global_maxsize The memory allocated.
#' @return A Seurat object containing the full dataset and a Seurat object containing only tumor cells.
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
  cnv_status = c("tumor", "normal"),
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

  # tumor_cells <- colnames(sc_seurat)[
  #   tolower(sc_seurat@meta.data$cnv_status) == cnv_status
  # ]
  # scRNA_tumor_dataset <- sc_seurat[, tumor_cells]

  scRNA_tumor_dataset = sc_seurat %>%
    subset(subset = tolower(.@meta.data$cnv_status) == cnv_status)

  result_list <- list(
    sc_seurat = process_sc_object(sc_seurat),
    sc_tumor_seurat = process_sc_object(scRNA_tumor_dataset)
  )

  return(result_list)
}


#' @description
#' Preprocess bulk expression data: Use the `IDconverter` package to convert Ensembles version IDs and TCGA version IDs to genes.
#' @param data raw bulk expression data
BulkPreProcess = function(data = NULL) {
  rownames(data) <- substr(rownames(data), 1, 15)
  options(
    IDConverter.datapath = system.file("extdata", package = "IDConverter")
  )
  rownames(data) <- IDConverter::convert_hm_genes(rownames(data))
  return(data)
}


#' @description
#' Preprocess mutational signature data
#' @param thresh Search for 6 types of mutational signatures. For each column, if the number of non-zero values is greater than `thresh`, it will be retained.
#' @param ms_search_pattern The matching pattern used to search for mutational signature columns.
MSPreProcess = function(
  ms_signature,
  thresh = 0.05,
  ms_search_pattern = "SBS|DBS|CN|CNV|SV|ID|INDEL"
) {
  library(dplyr)

  keep_colnames <- purrr::map_lgl(names(ms_signature), function(col_name) {
    col <- ms_signature[[col_name]]

    if (is.character(col)) {
      cli::cli_alert_info("{crayon::green(col_name)} kept")
      return(TRUE) # keep annotational columns
    } else if (is.numeric(col)) {
      ifelse(
        !grepl(ms_search_pattern, col_name),
        {
          cli::cli_alert_info(
            "{crayon::red(col_name)} not kept, because it represents non-mutational signature profile"
          )
          return(FALSE) # remove other unrelated data columns
        },
        ifelse(
          # Retain cols where the proportion of non-zero data in the ms exceeds the threshold.
          mean(col != 0) >= thresh, # equal to: sum(col != 0)/length(col) >= thresh
          {
            cli::cli_alert_info("{crayon::green(col_name)} kept")
            return(TRUE)
          },
          {
            cli::cli_alert_info(
              "{crayon::red(col_name)} not kept, because the valid data in the column < {crayon::bold(thresh)}"
            )
            return(FALSE)
          }
        )
      )
    } else {
      FALSE
    }
  })

  processed_ms_signature <- ms_signature[, keep_colnames, drop = FALSE]

  # final check for data availability
  if (sum(grepl(ms_search_pattern, colnames(processed_ms_signature))) == 0) {
    cli::cli_alert_warning(crayon::yellow(glue::glue(
      "{crayon::bold('All data columns were filtered out, possible reasons:')}",
      "1. Threshold too high (current: {crayon::black(thresh)})",
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

  return(list(
    phenotype = processed_ms_signature %>%
      dplyr::filter(.[[sample_colname]] %in% cm_samples) %>%
      dplyr::arrange(factor(.[[sample_colname]], levels = cm_samples)) %>%
      {
        stats::setNames(.$ms_status, .[[sample_colname]])
      },
    matched_TCGA_exp_count = TCGA_exp_count[, cm_samples, drop = FALSE],
    ms_select = names(ms_signature)[[col_id]]
  ))
}

# ------------------- other function ----------------------

# IdentifyDataColumn <- function(data, thresh = 0.5) {
#   pattern_based <- sapply(data, function(x) {
#     if (is.character(x)) {
#       numeric_pat <- grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", x)
#       mean(numeric_pat, na.rm = TRUE) > 0.5
#     } else {
#       TRUE
#     }
#   })
#   return(pattern_based)
# }
