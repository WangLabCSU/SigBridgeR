# * ---- Preprocess bulk expression data ----

#' @title Bulk RNA-seq Data Preprocessing and Quality Control Function
#'
#' @description
#' This function performs comprehensive preprocessing and quality control analysis
#' for bulk RNA-seq data, including data validation, filtering, batch effect detection,
#' principal component analysis, and visualization.
#'
#' @details
#' The function performs the following operations:
#' \enumerate{
#'   \item Data validation and format conversion
#'   \item Basic statistics calculation (missing values, read depth, gene detection)
#'   \item Optional detailed quality checks including:
#'     \itemize{
#'       \item Sample correlation analysis
#'       \item Principal Component Analysis (PCA)
#'       \item Outlier detection using Mahalanobis distance
#'       \item Batch effect detection using ANOVA
#'     }
#'   \item Data filtering based on count thresholds and quality metrics
#'   \item Optional gene symbol conversion
#'   \item Visualization generation (PCA plots)
#' }
#'
#' @param data Expression matrix with genes as rows and samples as columns, or a list containing count_matrix and sample_info
#' @param sample_info Sample information data frame (optional), ignored if data is a list. A qualified `sample_info` should contain both `sample` and `condition` columns (case-sensitive), and there are no specific requirements for the data type stored in the `condition` column. `batch` column is optional, which is used for batch effect detection.
#' @param gene_symbol_conversion Whether to convert Ensembles version IDs and TCGA version IDs to genes with IDConverter, default: `FALSE`
#' @param check Whether to perform detailed quality checks, default: `TRUE`
#' @param min_count_threshold Minimum count threshold for gene filtering. Only values greater than this threshold are considered to represent valid gene expression. Default: `10L`
#' @param min_gene_expressed Minimum number of samples a gene must be expressed in, default: `3L`
#' @param min_total_reads Minimum total reads per sample, default: `1e6L`
#' @param min_genes_detected Minimum number of genes detected per sample, default: `10000`
#' @param min_correlation Minimum correlation threshold between samples, default: `0.8`
#' @param n_top_genes Number of top variable genes for PCA analysis, default: `500`
#' @param show_plot_results Whether to generate visualization plots, default: `TRUE`
#' @param ... Additional arguments. Currently supports:
#'    - `verbose`: Logical indicating whether to print progress messages. Defaults to `TRUE`.
#'    - `seed`: For reproducibility, default is `123L`
#'    - `method`: Method for duplicates aggregation, see [AggregateDups]
#'
#' @section Quality Metrics:
#' The function calculates and reports several quality metrics:
#' \describe{
#'   \item{Data Integrity}{Number of missing values}
#'   \item{Gene Count}{Total number of genes after filtering}
#'   \item{Sample Read Depth}{Total reads per sample}
#'   \item{Gene Detection Rate}{Number of genes detected per sample}
#'   \item{Sample Correlation}{Pearson correlation between samples}
#'   \item{PCA Variance}{Variance explained by first two principal components}
#'   \item{Batch Effects}{Proportion of genes significantly affected by batch}
#' }
#'
#'
#' @section Sample Information Format:
#' The sample_info data frame should contain:
#' \describe{
#'   \item{sample}{Character vector of unique sample identifiers}
#'   \item{condition}{Character vector of experimental conditions}
#'   \item{batch}{Optional character vector of batch identifiers}
#' }
#'
#' @section Filtering Criteria:
#' Genes are retained if:
#' \itemize{
#'   \item They have counts >= min_count_threshold in >= min_gene_expressed samples
#' }
#'
#' Samples are retained if:
#' \itemize{
#'   \item Total reads >= min_total_reads
#'   \item Detected genes >= min_genes_detected
#'   \item Mean correlation >= min_correlation (if check=TRUE)
#' }
#'
#'
#' @seealso
#' \code{\link[edgeR]{cpm}} for counts per million calculation,
#' \code{\link[stats]{prcomp}} for PCA analysis,
#' \code{\link[stats]{cor}} for correlation analysis,
#' \code{\link[SigBridgeR]{SymbolConvert}} for gene symbol conversion
#'
#' @return Filtered count matrix
#' @export
#'
BulkPreProcess <- function(
  data,
  sample_info = NULL,
  gene_symbol_conversion = FALSE,
  check = TRUE,
  min_count_threshold = 10L,
  min_gene_expressed = 3L,
  min_total_reads = 1e6L,
  min_genes_detected = 10000L,
  min_correlation = 0.8,
  n_top_genes = 500L,
  show_plot_results = TRUE,
  ...
) {
  purrr::walk(
    list(
      min_count_threshold,
      min_gene_expressed,
      min_genes_detected,
      min_total_reads,
      min_correlation,
      n_top_genes
    ),
    ~ chk::chk_numeric
  )
  purrr::walk(
    list(gene_symbol_conversion, check, show_plot_results),
    ~ chk::chk_flag
  )

  # dots arguments
  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption("verbose")
  seed <- dots$seed %||% SigBridgeRUtils::getFuncOption("seed")
  method <- dots$method # Duplicate handling method

  set.seed(seed)

  # * Handle input data format
  if (is.list(data)) {
    if ("count_matrix" %chin% names(data)) {
      counts_matrix <- as.matrix(data$count_matrix)
      if ("sample_info" %chin% names(data) && is.null(sample_info)) {
        sample_info <- data$sample_info
      }
    } else {
      cli::cli_abort(c(
        "x" = "If data is a list, it must contain 'count_matrix' element"
      ))
    }
  } else {
    counts_matrix <- as.matrix(data)
  }

  n_genes <- nrow(counts_matrix)
  n_samples <- ncol(counts_matrix)

  if (verbose) {
    ts_cli$cli_alert_success(
      "Data loaded: {.val {n_genes}} genes * {.val {n_samples}} samples"
    )
  }

  if (gene_symbol_conversion) {
    counts_matrix <- SymbolConvert(counts_matrix)
    if (verbose) {
      ts_cli$cli_alert_success("Gene symbol conversion done")
    }
  }

  counts_matrix <- BulkCheck(
    counts_matrix = counts_matrix,
    n_genes = n_genes,
    min_genes_detected = min_genes_detected,
    min_count_threshold = min_count_threshold,
    method = method,
    verbose = verbose
  )

  if (is.null(sample_info)) {
    # Create default sample information
    if (verbose) {
      ts_cli$cli_alert_info(
        "No sample info provided, using default settings."
      )
    }
    sample_info <- data.frame(
      sample = colnames(counts_matrix) %||%
        glue::glue("Sample_{seq_len(n_samples)}"),
      condition = if (rlang::is_installed("cheapr")) {
        cheapr::rep_("unknown", n_samples)
      } else {
        rep("unknown", n_samples)
      },
      stringsAsFactors = FALSE
    )
  } else {
    # Validate matching
    if (nrow(sample_info) != n_samples) {
      cli::cli_abort(c(
        "x" = "Number of rows in `sample_info` does not match number of columns in count matrix"
      ))
    }
  }

  # * Basic Statistics
  missing_data <- sum(is.na(counts_matrix))
  # Data integrity
  if (missing_data != 0) {
    cli::cli_warn(
      "Data integrity: {.val {missing_data}} missing values detected"
    )
  }

  total_reads_per_sample <- SigBridgeRUtils::colSums3(
    counts_matrix,
    na.rm = TRUE
  )
  genes_detected_per_sample <- SigBridgeRUtils::colSums3(
    counts_matrix > 0,
    na.rm = TRUE
  )
  genes_with_reads <- SigBridgeRUtils::rowSums3(counts_matrix > 0, na.rm = TRUE)
  genes_expressed_multiple <- sum(genes_with_reads >= min_gene_expressed)

  # Detailed Quality Checks (Conditional Execution)
  if (check) {
    if (verbose) {
      ts_cli$cli_alert_info(
        "Starting detailed quality checks..."
      )
    }

    # Calculate CPM (compute once, reuse)
    cpm_values <- edgeR::cpm(counts_matrix, log = TRUE, prior.count = 1)

    # Sample correlation
    sample_cor <- if (rlang::is_installed("WGCNA")) {
      WGCNA::cor(cpm_values, method = "pearson", use = "pairwise.complete.obs")
    } else {
      stats::cor(cpm_values, method = "pearson", use = "pairwise.complete.obs")
    }

    # Quick identification of low correlations
    cor_lower_tri <- sample_cor[lower.tri(sample_cor)]
    min_cor <- min(cor_lower_tri, na.rm = TRUE)
    mean_cor <- mean(cor_lower_tri, na.rm = TRUE)
    low_cor_count <- sum(cor_lower_tri < min_correlation, na.rm = TRUE)

    if (verbose) {
      if (min_cor >= min_correlation) {
        cli::cli_alert_success(
          "Sample correlation: {cli::col_green('Good')} (minimum = {.val {round(min_cor, 3)}})"
        )
      } else {
        cli::cli_warn(
          "Sample correlation: {cli::col_red('Low')} (minimum = {.val {round(min_cor, 3)}})"
        )
      }
      if (low_cor_count > 0) {
        cli::cli_warn(
          "Found {.val {low_cor_count}} sample pairs with correlation < {.val min_correlation = {min_correlation}}"
        )
      }
    }

    # * PCA Analysis
    # Select highly variable genes
    n_genes_for_pca <- min(n_top_genes, nrow(cpm_values))
    gene_vars <- SigBridgeRUtils::rowVars3(cpm_values, na.rm = TRUE)
    top_var_idx <- order(gene_vars, decreasing = TRUE)[seq_len(
      n_genes_for_pca
    )]

    # PCA computation
    pca_data <- Matrix::t(cpm_values[top_var_idx, , drop = FALSE])
    pca_result <- stats::prcomp(pca_data, scale. = TRUE, center = TRUE)
    pca_summary <- summary(pca_result)

    # Outlier detection (use only first 2 PCs)
    if (ncol(pca_result$x) >= 2) {
      pc12 <- pca_result$x[, c(1, 2)]
      mahal_dist <- stats::mahalanobis(
        pc12,
        center = SigBridgeRUtils::colMeans3(pc12),
        cov = stats::cov(pc12)
      )
      outlier_threshold <- stats::qchisq(0.95, df = 2)
      outliers <- if (rlang::is_installed("cheapr")) {
        cheapr::which_(mahal_dist > outlier_threshold)
      } else {
        which(mahal_dist > outlier_threshold)
      }
    } else {
      outliers <- integer(0)
    }

    if (verbose) {
      var_pc1 <- round(pca_summary$importance[2, 1] * 100, 2)
      var_pc2 <- round(pca_summary$importance[2, 2] * 100, 2)
      ts_cli$cli_alert_info(
        "PCA completed: PC1({.val {var_pc1}}%) PC2({.val {var_pc2}}%), {.val {length(outliers)}} outlier samples"
      )
      if (length(outliers) > 0) {
        outlier_names <- sample_info$sample[outliers]
        cli::cli_warn(
          "Detected {.val {length(outliers)}} outlier sample(s) : {.emph {outlier_names}}"
        )
      }
    }

    # * Batch effect detection (ANOVA)
    if ("batch" %chin% colnames(sample_info)) {
      batch_levels <- unique(sample_info$batch)
      if (length(batch_levels) > 1) {
        if (verbose) {
          ts_cli$cli_alert_info(
            "Detecting batch effects..."
          )
        }

        nrow_cpm <- nrow(cpm_values)

        # Fast batch effect detection (random sampling of genes)
        batch_test_idx <- sample(nrow_cpm, min(5000, nrow_cpm))

        perform_anova_test <- function(gene_expression, batch_info) {
          if (all(is.na(gene_expression))) {
            return(1)
          }

          aov_result <- stats::aov(
            gene_expression ~ batch_info
          )
          summary(aov_result)[[1]]["Pr(>F)"][1, 1]
        }

        batch_pvalues <- cpm_values[batch_test_idx, , drop = FALSE] %>%
          asplit(1) %>% # by row
          purrr::map_dbl(~ perform_anova_test(.x, sample_info$batch))

        batch_sig_prop <- mean(batch_pvalues < 0.05, na.rm = TRUE)

        if (verbose) {
          ts_cli$cli_alert_success(
            "Batch effect detection completed: {.val {round(batch_sig_prop * 100, 2)}}% genes affected"
          )
          if (batch_sig_prop > 0.3) {
            cli::cli_warn(
              "Strong batch effects detected, consider batch correction."
            )
          }
        }
      }
    }

    # Generate key plots
    if (show_plot_results) {
      if (verbose) {
        ts_cli$cli_alert_info(
          "Generating visualization plots..."
        )
      }
      rlang::check_installed(c("ggplot2", "ggforce"))

      # PCA graphics::plot
      if (exists("pca_result")) {
        pca_df <- data.frame(
          PC1 = pca_result$x[, 1],
          PC2 = pca_result$x[, 2],
          sample = sample_info$sample,
          condition = sample_info$condition,
          stringsAsFactors = FALSE
        )

        if ("batch" %chin% colnames(sample_info)) {
          pca_df$batch <- sample_info$batch
        }

        var_labels <- pca_summary$importance[2, c(1, 2)] * 100

        p_pca <- ggplot2::ggplot(
          pca_df,
          ggplot2::aes(x = `PC1`, y = `PC2`, color = `condition`)
        ) +
          ggplot2::geom_point(size = 3, alpha = 0.8) +
          ggplot2::labs(
            title = "Principal Component Analysis (PCA)",
            x = paste0("PC1 (", round(var_labels[1], 2), "%)"),
            y = paste0("PC2 (", round(var_labels[2], 2), "%)")
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = "right") +
          ggforce::geom_mark_ellipse(
            ggplot2::aes(fill = condition, group = condition),
            alpha = 0.1,
            expand = ggplot2::unit(3, "mm"),
            show.legend = FALSE
          )

        if ("batch" %chin% colnames(pca_df)) {
          p_pca <- p_pca +
            ggplot2::aes(shape = `batch`) +
            ggplot2::scale_shape_manual(
              values = c(
                16,
                17,
                18,
                19,
                20
              )[seq_len(length(unique(pca_df$batch)))]
            )
        }

        methods::show(p_pca)
      }
    }

    rm(cpm_values, pca_data, pca_result)
    if (exists("pc12")) {
      rm(pc12)
    }
    gc(verbose = FALSE)
  }

  # * Data Filtering
  # Gene filtering
  genes_pass_count <- SigBridgeRUtils::rowSums3(
    counts_matrix >= min_count_threshold
  ) >=
    min_gene_expressed

  # Sample filtering
  samples_pass_reads <- total_reads_per_sample >= min_total_reads # will be used later
  samples_pass_genes <- genes_detected_per_sample >= min_genes_detected # will be used later
  samples_to_keep <- samples_pass_reads & samples_pass_genes

  # Correlation filtering (if correlation was computed)
  if (check && exists("sample_cor")) {
    sample_mean_cor <- SigBridgeRUtils::colMeans3(sample_cor, na.rm = TRUE)
    samples_to_keep <- samples_to_keep & sample_mean_cor >= min_correlation
  }

  filtered_counts <- counts_matrix[
    genes_pass_count,
    samples_to_keep,
    drop = FALSE
  ]

  if (verbose) {
    n_genes_filtered <- sum(genes_pass_count)
    n_samples_filtered <- sum(samples_to_keep)

    genes_removed <- n_genes - n_genes_filtered
    samples_removed <- n_samples - n_samples_filtered
    ts_cli$cli_alert_success("Data filtering completed:")
    cli::cli_alert_info(
      "  Genes: {.val {n_genes}} -> {.val {n_genes_filtered}} (removed {.val {genes_removed}})"
    )
    cli::cli_alert_info(
      "  Samples: {.val {n_samples}} -> {.val {n_samples_filtered}} (removed {.val {samples_removed}})"
    )

    # * Quality Report Generation
    if (!all(samples_pass_reads[samples_to_keep])) {
      n_low_reads <- sum(!samples_pass_reads[samples_to_keep])
      cli::cli_warn(
        "Sample read depth: {.val {n_low_reads}} samples below threshold"
      )
    }

    if (!all(samples_pass_genes[samples_to_keep])) {
      n_low_genes <- sum(!samples_pass_genes[samples_to_keep])
      cli::cli_warn(
        "Gene detection: {.val {n_low_genes}} samples below threshold"
      )
    }

    ts_cli$cli_alert_success(
      cli::col_green("BulkPreProcess completed")
    )
  }

  filtered_counts
}

#' @keywords internal
BulkCheck <- function(
  counts_matrix,
  n_genes,
  min_genes_detected,
  min_count_threshold,
  method,
  verbose = TRUE
) {
  if (any(counts_matrix < 0, na.rm = TRUE)) {
    cli::cli_abort(c(
      "x" = "Expression matrix cannot contain negative values"
    ))
  }
  if (any(counts_matrix != floor(counts_matrix))) {
    cli::cli_warn("Expression matrix contains non-integre value")
  }
  if (min_genes_detected > n_genes) {
    cli::cli_abort(c(
      "x" = "{.arg min_genes_detected} must be less than the number of genes in the data",
      ">" = "Current number of genes: {.val {n_genes}}"
    ))
  }
  if (min_count_threshold > max(counts_matrix)) {
    cli::cli_abort(c(
      "x" = "{.arg min_count_threshold} must be less than the maximum count in the data",
      ">" = "Current maximum count: {.val {max(counts_matrix)}}"
    ))
  }

  # * Handle duplicated genes and samples
  if (any(duplicated(rownames(counts_matrix)))) {
    cli::cli_alert_info("Aggregate Duplicated genes in rownames")
    counts_matrix <- AggregateDupRows(
      counts_matrix,
      method = method,
      verbose = FALSE
    )
  }
  if (any(duplicated(colnames(counts_matrix)))) {
    cli::cli_alert_info("Aggregate Duplicated samples in colnames")
    counts_matrix <- AggregateDupCols(
      counts_matrix,
      method = method,
      verbose = FALSE
    )
  }

  # * Handle data type
  if (!is.numeric(counts_matrix)) {
    cli::cli_warn(
      "Expression matrix must be numeric, converted to numeric now."
    )
    # Attempt to convert to numeric
    counts_matrix <- as.numeric(counts_matrix)
    # Check if conversion introduced NAs
    new_nas <- sum(is.na(counts_matrix))
    if (new_nas > 0) {
      cli::cli_warn(
        "Numeric conversion introduced {.val {new_nas}} NA values (will be imputed)"
      )
    }
  }

  counts_matrix
}
