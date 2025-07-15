# * ---- scRNA-seq preprocessing ----

#' @title Single-Cell RNA-seq Preprocessing Pipeline
#'
#' @description
#' A generic function for standardized preprocessing of single-cell RNA-seq data
#' from multiple sources. Handles data.frame/matrix, AnnData, and Seurat inputs
#' with tumor cell filtering. Implements a complete analysis pipeline
#' from raw data to clustered embeddings.
#'
#' @name SCPreProcess
#' @usage
#' SCPreProcess(sc, ...)
#'
#' @param sc Input data, one of:
#'        - `data.frame/matrix`: Raw count matrix (features x cells)
#'        - `AnnDataR6`: Python AnnData object via reticulate
#'        - `Seurat`: Preprocessed Seurat object
#' @param ... Method-specific arguments (see below)
#'
#' @return A Seurat object containing:
#' \itemize{
#'   \item Normalized and scaled expression data
#'   \item Variable features
#'   \item PCA/tSNE/UMAP reductions
#'   \item Cluster identities
#'   \item When tumor cells filtered: original dimensions in `@misc$raw_dim`
#'   \item Final dimensions in `@misc$self_dim`
#' }
NULL


#' @rdname SCPreProcess
#' @export
SCPreProcess <- function(sc, ...) {
    UseMethod("SCPreProcess")
}

#' @rdname SCPreProcess
#' @export
SCPreProcess.default <- function(sc, ...) {
    stop("Unknown Input type")
}

#' @rdname SCPreProcess
#' @param column2only_tumor Metadata column for tumor cell filtering (regex patterns:
#'        "[Tt]umo.?r", "[Cc]ancer", "[Mm]alignant", "[Nn]eoplasm")
#' @param project Project name for Seurat object
#' @param min_cells Minimum cells per gene to retain
#' @param min_features Minimum features per cell to retain
#' @param normalization_method Normalization method ("LogNormalize", "CLR", or "RC")
#' @param scale_factor Scaling factor for normalization
#' @param selection_method Variable feature selection method ("vst", "mvp", or "disp")
#' @param resolution Cluster resolution (higher for more clusters)
#' @param dims PCA dimensions to use
#' @param verbose Print progress messages
#' @param future_global_maxsize Memory limit for parallelization (bytes)
#' @export
#'
SCPreProcess.data.frame <- function(
    sc,
    column2only_tumor = NULL,
    project = "Scissor_Single_Cell",
    min_cells = 400,
    min_features = 0,
    normalization_method = "LogNormalize",
    scale_factor = 10000,
    selection_method = "vst",
    resolution = 0.6,
    dims = 1:10,
    verbose = TRUE,
    future_global_maxsize = 6 * 1024^3,
    ...
) {
    options(future.globals.maxSize = future_global_maxsize)

    # sc is a count matrix
    cli::cli_alert_info("Start from count matrix")

    sc_seurat <- SeuratObject::CreateSeuratObject(
        counts = sc,
        project = project,
        min.cells = min_cells,
        min.features = min_features
    ) %>%
        ProcessSeuratObject(
            normalization_method = normalization_method,
            scale_factor = scale_factor,
            selection_method = selection_method,
            verbose = verbose
        )

    # Add metadata
    if ("obs" %in% names(sc)) {
        sc_seurat <- sc_seurat %>% Seurat::AddMetaData(sc$obs)
    }

    sc_seurat <- ClusterAndReduce(
        sc_seurat,
        dims = dims,
        resolution = resolution,
        verbose = verbose
    )

    FilterTumorCell(
        sc_seurat,
        column2only_tumor = column2only_tumor,
        verbose = verbose
    )
}

#' @rdname SCPreProcess
#' @export
SCPreProcess.AnnDataR6 <- function(
    sc,
    column2only_tumor = NULL,
    project = "Scissor_Single_Cell",
    min_cells = 400,
    min_features = 0,
    normalization_method = "LogNormalize",
    scale_factor = 10000,
    selection_method = "vst",
    resolution = 0.6,
    dims = 1:10,
    verbose = TRUE,
    future_global_maxsize = 6 * 1024^3,
    ...
) {
    options(future.globals.maxSize = future_global_maxsize)

    if (is.null(sc$X)) {
        stop("Input must contain $X matrix")
    }
    cli::cli_alert_info("Start from anndata object")

    sc_matrix <- if (inherits(sc$X, "sparseMatrix")) {
        if (!requireNamespace("Matrix", quietly = TRUE)) {
            stop("Matrix package required for sparse matrix support")
        }
        Matrix::t(sc$X)
    } else {
        t(sc$X)
    }

    sc_seurat <- SeuratObject::CreateSeuratObject(
        counts = sc_matrix,
        project = project,
        min.cells = min_cells,
        min.features = min_features
    ) %>%
        ProcessSeuratObject(
            normalization_method = normalization_method,
            scale_factor = scale_factor,
            selection_method = selection_method,
            verbose = verbose
        )

    # Add metadata
    if ("obs" %in% names(sc)) {
        sc_seurat <- sc_seurat %>% Seurat::AddMetaData(sc$obs)
    }

    sc_seurat <- ClusterAndReduce(
        sc_seurat,
        dims = dims,
        resolution = resolution,
        verbose = verbose
    )

    FilterTumorCell(
        sc_seurat,
        column2only_tumor = column2only_tumor,
        verbose = verbose
    )
}

#' @rdname SCPreProcess
#' @export
#'
SCPreProcess.Seurat <- function(
    sc,
    column2only_tumor = NULL,
    verbose = TRUE,
    future_global_maxsize = 6 * 1024^3,
    ...
) {
    options(future.globals.maxSize = future_global_maxsize)

    FilterTumorCell(
        obj = sc,
        column2only_tumor = column2only_tumor,
        verbose = verbose
    )
}

#' @title Process a Seurat object (internal)
#'
#' @description
#' Normalize, find variable features, scale, and run PCA
#'
#' @param obj Seurat object
#' @param normalization_method Normalization method ("LogNormalize", "CLR", or "RC")
#' @param scale_factor Scaling factor for normalization
#' @param selection_method Variable feature selection method ("vst", "mvp", or "disp")
#' @param verbose Print progress messages
#' @return Seurat object
#'
#' @keywords internal
#'
ProcessSeuratObject <- function(
    obj,
    normalization_method = "LogNormalize",
    scale_factor = 10000,
    selection_method = "vst",
    verbose = TRUE
) {
    obj %>%
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
        )
}

#' Cluster and reduce dimensions (internal)
#'
#' @description
#' FindNeighbors, FindClusters, RunTSNE, RunUMAP
#'
#'
#' @keywords internal
#'
ClusterAndReduce <- function(
    obj,
    dims = 1:10,
    resolution = 0.6,
    verbose = TRUE
) {
    n_pcs <- ncol(obj@reductions$pca)
    if (is.null(dims) || max(dims) > n_pcs) {
        dims <- 1:min(max(dims), n_pcs)
        cli::cli_alert_warning(crayon::yellow(
            "The input dimension is greater than the dimension in PCA. It is now set to the maximum dimension in PCA."
        ))
    }

    obj %>%
        Seurat::FindNeighbors(dims = dims, verbose = verbose) %>%
        Seurat::FindClusters(
            resolution = resolution,
            verbose = verbose
        ) %>%
        Seurat::RunTSNE(dims = dims) %>%
        Seurat::RunUMAP(dims = dims, verbose = verbose)
}

#' Filter tumor cells (internal)
#'
#' @description
#' Filter tumor cells from Seurat object.
#'
#' @param obj Seurat object with a column to filter out tumor cells.
#' @param name2only_tumor Name of the column to filter out tumor cells.
#' @param verbose Logical. Whether to print messages.
#'
#' @keywords internal
#'
#'
FilterTumorCell <- function(
    obj,
    column2only_tumor = NULL,
    verbose = TRUE
) {
    obj = AddMisc(obj, self_dim = dim(obj), cover = TRUE)

    if (!is.null(column2only_tumor)) {
        ifelse(
            !column2only_tumor %in% colnames(obj@meta.data),
            {
                cli::cli_alert_danger(crayon::red(
                    "Column '{column2only_tumor}' not found, skip tumor cell filtering"
                ))
                return(obj)
            },
            {
                labels <- obj[[column2only_tumor]][[1]]
                tumor_cells <- grepl(
                    "^[Tt]umo.?r|[Cc]ancer[Mm]alignant|[Nn]eoplasm",
                    labels
                )

                tumor_seurat <- obj[, tumor_cells] %>%
                    AddMisc(
                        raw_dim = dim(obj),
                        self_dim = dim(.),
                        column2only_tumor = column2only_tumor,
                        cover = TRUE
                    )

                return(list(tumor_seurat = tumor_seurat, raw_seurat = obj))
            }
        )
    } else {
        return(obj)
    }
}


# * ---- Preprocess bulk expression data ----

#' @title Preprocess bulk expression data
#'
#' @description
#' Preprocess bulk expression data: convert Ensembles version IDs and TCGA version IDs to genes.
#' @param data raw bulk expression data
#'
#' @export
#'
BulkPreProcess = function(data) {
    #   rownames(data) <- substr(rownames(data), 1, 15)
    options(
        IDConverter.datapath = system.file("extdata", package = "IDConverter")
    )
    rownames(data) <- IDConverter::convert_hm_genes(rownames(data))
    return(data)
}


# * ---- Preprocess Mutational Signature Data ----

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
#'
MSPreProcess <- function(
    ms_signature,
    filter_tumor_type = "all",
    col_thresh = 0.05,
    accuracy_thresh = 0,
    ms_search_pattern = "SBS|DBS|CN|CNV|SV|ID|INDEL"
) {
    library(dplyr)

    # *filter columns under threshold
    keep_colnames <- lapply(X = names(ms_signature), FUN = function(col_name) {
        col <- ms_signature[[col_name]]

        if (is.character(col)) {
            cli::cli_alert_info("{crayon::green(col_name)} kept")
            return(col_name) # keep annotational columns
        } else if (is.numeric(col)) {
            if (!grepl(tolower(ms_search_pattern), tolower(col_name))) {
                cli::cli_alert_info(
                    "{crayon::red(col_name)} not kept, because it does not match the pattern"
                )
                return(NULL)
            }

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
    ms_signature <- ms_signature %>%
        dplyr::select(unlist(keep_colnames))

    # *filter accuracy
    accuracy_column <- grep(
        "[aA][cC]{2}.*",
        colnames(ms_signature),
        value = TRUE
    )
    if (length(accuracy_column) == 0) {
        cli::cli_alert_warning(crayon::yellow(
            "No accuracy column found in the data"
        ))
    } else {
        ms_signature <- ms_signature %>%
            dplyr::filter(
                .[[accuracy_column]] >= accuracy_thresh,
                grepl(
                    pattern = glue::glue(
                        tolower(filter_tumor_type),
                        .sep = "|"
                    ),
                    x = tolower(.[[tumor_type_col]])
                )
            )
    }
    # *filter cancer type
    tumor_type_col <- grep(
        "[tT]umo.?r|[Cc]ancer",
        colnames(ms_signature),
        value = TRUE
    )
    if (length(tumor_type_col) == 0) {
        cli::cli_alert_warning(crayon::yellow(
            "No tumor type column found in the data"
        ))
    } else {
        all_tumor_types = unique(ms_signature[[tumor_type_col]])
        cli::cli_alert_info(glue::glue(
            crayon::bold("Tumor types in data:"),
            all_tumor_types
        ))
        cli::cli_alert_info(c(
            crayon::bold("Filtering by:"),
            " {filter_tumor_type}"
        ))

        if (
            tolower(filter_tumor_type) %in% c("all", "all tumor", "all types")
        ) {
            filter_tumor_type = all_tumor_types
        }

        ms_signature <- ms_signature %>%
            dplyr::filter(
                grepl(
                    pattern = glue::glue(
                        tolower(filter_tumor_type),
                        .sep = "|"
                    ),
                    x = tolower(.[[tumor_type_col]])
                )
            )
    }

    # final check for data availability
    if (sum(grepl(ms_search_pattern, colnames(ms_signature))) == 0) {
        cli::cli_alert_warning(crayon::yellow(glue::glue(
            "{crayon::bold('All data columns were filtered out, possible reasons:')}",
            "1. Threshold too high (current: {crayon::black(col_thresh)})",
            "2. No matching data columns existed (search pattern: {crayon::black(ms_search_pattern)})",
            "3. All values were zero",
            .sep = "\n"
        )))
    }

    return(ms_signature)
}


#' @title Match Samples Across TCGA Datasets
#'
#' @description
#' Identifies common samples between mutational signature data and TCGA expression
#' count matrix by sample ID truncation, then returns matched subsets with binary
#' phenotype encoding. Designed to resolve inconsistencies in TCGA sample identifiers
#' across different genomic datasets.
#'
#' @param ms_signature A data frame containing mutational signature data. Must include:
#'    - A sample ID column (auto-detected via case-insensitive "sample" pattern)
#'    - A numeric column specified by `col_id` for binarization.
#' @param bulk_data A matrix or data frame of expression counts with samples
#'    as columns. Column names should contain TCGA sample IDs.
#' @param col_id Integer, string or vector. Column index/name(s) in `ms_signature` to use for:
#'    - Binarization (values > `ms_status_thresh` become 1, others 0)
#'    - Reporting in output as `ms_select`
#'    - When multiple mutational signatures are specified, the `ms_status` will only be recorded as `1` if each mutational signature > `ms_status_thresh`.
#' 
#' @param ms_status_thresh Numeric. Threshold for converting the `col_id` column to
#'    binary status (default: 0L). Values above threshold become 1.
#'
#' @return A list with three components:
#' \describe{
#'   \item{phenotype}{Named integer vector (names = sample IDs) of binarized values}
#'   \item{matched_bulk_data}{Expression matrix subset to common samples}
#'   \item{ms_select}{Character indicating which signature was used (column name)}
#' }
#'
#' @section Processing Details:
#' 1. **Sample ID Harmonization**:
#'    - Truncates IDs to the minimum length across both datasets (max 15 chars)
#'    - Case-sensitive exact matching after truncation
#' 2. **Phenotype Conversion**:
#'    - Converts `col_id` values to binary (0/1) based on threshold
#' 3. **Validation**:
#'    - Stops if no common samples found
#'    - Preserves sample order in output
#'
#' @examples
#' \dontrun{
#' # Mock data
#' sig_data <- data.frame(
#'   sample_id = c("TCGA-AB-1234-01", "TCGA-CD-5678-01"),
#'   APOBEC_score = c(0.8, 0.2)
#' )
#' exp_matrix <- matrix(rnorm(20), ncol = 4,
#'   dimnames = list(NULL, c("TCGA-AB-1234", "TCGA-EF-9012", "TCGA-CD-5678", "TCGA-GH-3456")))
#'
#' # Basic usage
#' result <- MatchSample(
#'   ms_signature = sig_data,
#'   bulk_data = exp_matrix,
#'   col_id = "APOBEC_score", # or 1L
#'   ms_status_thresh = 0.5
#' )
#'
#' # Output components
#' str(result)
#' # $ phenotype       : Named int [1:2] 1 0
#' # $ matched_TCGA_exp: num [1:5, 1:2] 0.481 -0.788 0.741 -0.311 -0.668 ...
#' # $ ms_select       : chr "APOBEC_score"
#' }
#'
#' @export
#' @importFrom dplyr mutate filter arrange
#' @importFrom cli cli_alert_info
#' @importFrom stats setNames
#'
MatchSample = function(
    ms_signature,
    bulk_data,
    col_id,
    ms_status_thresh = 0L
) {
    require(dplyr)
    if (
        !all(col_id %in% colnames(ms_signature)) &&
            !all(col_id %in% seq_along(ms_signature))
    ) {
        cli::cli_abort(
            c(
                "x" = "{.var col_id}={.val {col_id}} not found in {.var ms_signature}",
                "i" = "Please specify a valid column name or index"
            ),
            class = "NotFoundError"
        )
    }
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
        min(nchar(colnames(bulk_data))),
        15
    )
    colnames(bulk_data) = substr(colnames(bulk_data), 1, trunc_length)
    cli::cli_alert_info("Sample match length: {trunc_length}")

    # convert interested column to binary variable
    processed_ms_signature <- ms_signature %>%
        dplyr::mutate(
            # Samples exceeding the threshold will be retained and set logical
            !!sample_colname := substr(.[[sample_colname]], 1, trunc_length),
            ms_status = as.integer(
                if (is.numeric(col_id)) {
                    rowSums(.[, col_id, drop = FALSE] > ms_status_thresh) ==
                        length(col_id)
                } else {
                    rowSums(dplyr::across(
                        dplyr::any_of(col_id),
                        ~ .x > ms_status_thresh
                    )) ==
                        length(col_id)
                }
            )
        )
    # find common sample names
    cm_samples = intersect(
        processed_ms_signature[[sample_colname]],
        colnames(bulk_data)
    )
    cli::cli_alert_info("Sample match: {length(cm_samples)} common samples")

    if (length(cm_samples) == 0) {
        cli::cli_abort(
            c("x" = "No common sample found in data"),
            class = "NullCommonSample"
        )
    }

    match_result = list(
        phenotype = processed_ms_signature %>%
            dplyr::filter(.[[sample_colname]] %in% cm_samples) %>%
            dplyr::arrange(factor(.[[sample_colname]], levels = cm_samples)) %>%
            {
                stats::setNames(.$ms_status, .[[sample_colname]])
            },
        matched_bulk_data = bulk_data[, cm_samples, drop = FALSE],
        ms_select = names(ms_signature)[[col_id]]
    )

    return(match_result)
}
