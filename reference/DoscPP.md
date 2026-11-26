# Perform scPP screening analysis

This function performs scPP screening on single-cell data using matched
bulk data and phenotype information. It supports binary, continuous, and
survival phenotype types.

## Usage

``` r
DoscPP(
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
)
```

## Arguments

- matched_bulk:

  Matrix or data frame of preprocessed bulk RNA-seq expression data
  (genes x samples). Column names must match names/IDs in phenotype.

- sc_data:

  A Seurat object containing scRNA-seq data to be screened.

- phenotype:

  Phenotype data, either:

  - Named vector, where names match columns of matched_bulk

  - Patient survival data frame, where row names match columns of
    matched_bulk, and column names are "time" and "status"

- label_type:

  Character specifying phenotype label type, e.g. "SBS1" or "time".

- phenotype_class:

  Type of phenotypic outcome. Must be consistent with input data:

  - "binary": Binary traits, e.g. case/control

  - "continuous": Continuous measurements

  - "survival": Survival information

- ref_group:

  Reference group or baseline for **binary** comparisons, e.g. "Normal"
  for Tumor/Normal studies and 0 for 0/1 case-control studies. (default:
  0)

- Log2FC_cutoff:

  Minimum log2 fold-change for binary markers (default: 0.585)

- estimate_cutoff:

  Effect size threshold for **continuous** traits (default: 0.2)

- probs:

  A numeric value indicating the quantile cutoff for cell
  classification. This parameter can also be a numeric vector, in which
  case an optimal threshold will be selected based on the AUC and
  enrichment score.(default: `0.2`)

- ...:

  Additional arguments. Currently supports:

  - `verbose`: Logical indicating whether to print progress messages.
    Defaults to `TRUE`.

  - `seed`: For reproducibility, default is `123L`

  - `assay`: The assay to use for single-cell data. Defaults to `RNA`.

## Value

A list containing:

- **scRNA_data**: Seurat object with added metadata:

  - **ScPP**: "Positive"/"Negative"/"Neutral" classification

- **gene_list**: List of genes used for screening

## Algorithm Steps

1.  Data Validation: Checks sample alignment between bulk and phenotype
    data

2.  Marker Selection: Identifies phenotype-associated genes from bulk
    data

3.  Single-cell Screening: Projects bulk markers onto single-cell data

4.  Cell Classification: Categorizes cells based on phenotype
    association

## Reference

WangX-Lab/ScPP \[Internet\]. \[cited 2025 Aug 31\]. Available from:
https://github.com/WangX-Lab/ScPP

## See also

Other screen_method:
[`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md),
[`DoLP_SGL()`](https://wanglabcsu.github.io/sigbridger/reference/DoLP_SGL.md),
[`DoPIPET()`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md),
[`DoPhenoMapR()`](https://wanglabcsu.github.io/sigbridger/reference/DoPhenoMapR.md),
[`DoSCIPAC()`](https://wanglabcsu.github.io/sigbridger/reference/DoSCIPAC.md),
[`DoSIDISH()`](https://wanglabcsu.github.io/sigbridger/reference/DoSIDISH.md),
[`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md),
[`DoTiRank()`](https://wanglabcsu.github.io/sigbridger/reference/DoTiRank.md),
[`DoscAB()`](https://wanglabcsu.github.io/sigbridger/reference/DoscAB.md),
[`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md)

Other scPP:
[`ValidatescPPParams()`](https://wanglabcsu.github.io/sigbridger/reference/ValidatescPPParams.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Binary phenotype analysis
res <- DoscPP(
  matched_bulk = bulk_data,
  sc_data = seurat_obj,
  phenotype = ms_data,
  label_type = "SBS1",
  phenotype_class = "Binary"
)

# Survival analysis
surv_res <- DoscPP(
  sc_data = seurat_obj,
  matched_bulk = bulk_data,
  phenotype = surv_df,
  label_type = "OS_status",
  phenotype_class = "Survival"
)
} # }
```
