# Convert Ensembles Version IDs & TCGA Version IDs to Genes in Bulk Expression Data

Preprocess bulk expression data: convert Ensembles version IDs and TCGA
version IDs to genes. NA values are replaced with `unknown_k` format (k
stands for the position of the NA value in the row).

## Usage

``` r
SymbolConvert(
  data,
  unknown_format = "unknown_{k}",
  genome_build = c("hg38", "hg19", "mm10", "mm9"),
  verbose = SigBridgeRUtils::getFuncOption("verbose"),
  ...
)
```

## Arguments

- data:

  bulk expression data (matrix or data.frame)

- unknown_format:

  A glue pattern containing `{k}` for replace the NA value during
  conversion. k must be wrapped in curly braces, stands for the position
  of the NA value in the row. Default: `"unknown_{k}"`.

- genome_build:

  Genome build of the data. Default: `"hg38"`.

- verbose:

  Whether to print verbose messages. Default: `TRUE`.

- ...:

  No usage
