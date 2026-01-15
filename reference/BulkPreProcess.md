# Bulk RNA-seq Data Preprocessing and Quality Control Function

This function performs comprehensive preprocessing and quality control
analysis for bulk RNA-seq data, including data validation, filtering,
batch effect detection, principal component analysis, and visualization.

## Usage

``` r
BulkPreProcess(
  data,
  sample_info = NULL,
  gene_symbol_conversion = FALSE,
  check = TRUE,
  min_count_threshold = 10L,
  min_gene_expressed = 3L,
  min_total_reads = 1000000L,
  min_genes_detected = 10000L,
  min_correlation = 0.8,
  n_top_genes = 500L,
  show_plot_results = TRUE,
  ...
)
```

## Arguments

- data:

  Expression matrix with genes as rows and samples as columns, or a list
  containing count_matrix and sample_info

- sample_info:

  Sample information data frame (optional), ignored if data is a list. A
  qualified `sample_info` should contain both `sample` and `condition`
  columns (case-sensitive), and there are no specific requirements for
  the data type stored in the `condition` column. `batch` column is
  optional, which is used for batch effect detection.

- gene_symbol_conversion:

  Whether to convert Ensembles version IDs and TCGA version IDs to genes
  with IDConverter, default: `FALSE`

- check:

  Whether to perform detailed quality checks, default: `TRUE`

- min_count_threshold:

  Minimum count threshold for gene filtering. Only values greater than
  this threshold are considered to represent valid gene expression.
  Default: `10L`

- min_gene_expressed:

  Minimum number of samples a gene must be expressed in, default: `3L`

- min_total_reads:

  Minimum total reads per sample, default: `1e6L`

- min_genes_detected:

  Minimum number of genes detected per sample, default: `10000`

- min_correlation:

  Minimum correlation threshold between samples, default: `0.8`

- n_top_genes:

  Number of top variable genes for PCA analysis, default: `500`

- show_plot_results:

  Whether to generate visualization plots, default: `TRUE`

- ...:

  Additional arguments. Currently supports:

  - `verbose`: Logical indicating whether to print progress messages.
    Defaults to `TRUE`.

  - `seed`: For reproducibility, default is `123L`

  - `method`: Method for duplicates aggregation, see
    [AggregateDups](https://wanglabcsu.github.io/sigbridger/reference/aggregate-dups.md)

## Value

Filtered count matrix

## Details

The function performs the following operations:

1.  Data validation and format conversion

2.  Basic statistics calculation (missing values, read depth, gene
    detection)

3.  Optional detailed quality checks including:

    - Sample correlation analysis

    - Principal Component Analysis (PCA)

    - Outlier detection using Mahalanobis distance

    - Batch effect detection using ANOVA

4.  Data filtering based on count thresholds and quality metrics

5.  Optional gene symbol conversion

6.  Visualization generation (PCA plots)

## Quality Metrics

The function calculates and reports several quality metrics:

- Data Integrity:

  Number of missing values

- Gene Count:

  Total number of genes after filtering

- Sample Read Depth:

  Total reads per sample

- Gene Detection Rate:

  Number of genes detected per sample

- Sample Correlation:

  Pearson correlation between samples

- PCA Variance:

  Variance explained by first two principal components

- Batch Effects:

  Proportion of genes significantly affected by batch

## Sample Information Format

The sample_info data frame should contain:

- sample:

  Character vector of unique sample identifiers

- condition:

  Character vector of experimental conditions

- batch:

  Optional character vector of batch identifiers

## Filtering Criteria

Genes are retained if:

- They have counts \>= min_count_threshold in \>= min_gene_expressed
  samples

Samples are retained if:

- Total reads \>= min_total_reads

- Detected genes \>= min_genes_detected

- Mean correlation \>= min_correlation (if check=TRUE)

## See also

[`cpm`](https://rdrr.io/pkg/edgeR/man/cpm.html) for counts per million
calculation, [`prcomp`](https://rdrr.io/r/stats/prcomp.html) for PCA
analysis, [`cor`](https://rdrr.io/r/stats/cor.html) for correlation
analysis,
[`SymbolConvert`](https://wanglabcsu.github.io/sigbridger/reference/SymbolConvert.md)
for gene symbol conversion
