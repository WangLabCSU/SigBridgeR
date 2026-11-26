# Convert Ensembles Version IDs & TCGA Version IDs to Genes in Bulk Expression Data

Convert Ensembles version IDs and TCGA version IDs to genes. NA values
(unmatched IDs) are replaced with the original unconverted IDs (Ensembl
version IDs or TCGA version IDs).

## Usage

``` r
SymbolConvert(
  data,
  genome_build = c("hg38", "hg19", "mm10", "mm9", "ce11", "T2T"),
  update_symbol = TRUE,
  update_symbol_method = c("auto", "scCustomize", "Seurat"),
  verbose = SigBridgeRUtils::getFuncOption("verbose"),
  unknown_format = lifecycle::deprecated(),
  ...
)
```

## Arguments

- data:

  RNA expression data with rownames as Ensembles version IDs or TCGA
  version IDs(matrix or data.frame)

- genome_build:

  Genome build of the data. Default: `"hg38"`.

- update_symbol:

  Whether to update gene symbol to latest version. Default: `TRUE`.

- update_symbol_method:

  Method used to update gene symbols, one of `"auto"`, `"scCustomize"`
  or `"Seurat"`. `"auto"` prefers `scCustomize` when it is installed,
  otherwise falls back to `Seurat`. Only used when
  `update_symbol = TRUE`. Default: `"auto"`.

- verbose:

  Whether to print verbose messages. Default: `TRUE`.

- unknown_format:

  **\[deprecated\]** This argument is deprecated and ignored. NA values
  are now replaced with the original unconverted IDs.

- ...:

  No usage

## Value

Data with rownames converted to gene symbols

## See also

Other input_preprocess:
[`DrawBulkPCA()`](https://wanglabcsu.github.io/sigbridger/reference/DrawBulkPCA.md),
[`PhenoMap()`](https://wanglabcsu.github.io/sigbridger/reference/PhenoMap.md),
[`PhenoPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/PhenoPreProcess.md),
[`SCPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcess.md)
