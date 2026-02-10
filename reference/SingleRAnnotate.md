# Annotate Single-Cell Data Using SingleR

Performs automated cell type annotation using the `SingleR` package.
Supports built-in references (e.g., Human Primary Cell Atlas) or
user-provided `SingleCellExperiment` objects. Integrates seamlessly with
`Seurat` objects by adding prediction results as metadata columns.

## Usage

``` r
SingleRAnnotate(
  sc,
  verbose = getFuncOption("verbose"),
  ensembl = FALSE,
  cell.ont = c("all", "nonna", "none"),
  legacy = FALSE,
  ref = c("HPCA", "custom"),
  labels = NULL,
  ...
)
```

## Arguments

- sc:

  A single-cell dataset, either:

  - A `Seurat` object (recommended), or

  - A `SingleCellExperiment` object.

- verbose:

  Logical. Whether to print progress messages. Default: inherits from
  package option `getOption("SigBridgeR.verbose")`.

- ensembl:

  Logical scalar indicating whether to convert row names to Ensembl IDs.
  Genes without a mapping to a non-duplicated Ensembl ID are discarded.

- cell.ont:

  String specifying whether Cell Ontology terms should be included in
  the colData. If `"nonna"`, all samples without a valid term are
  discarded; if `"all"`, all samples are returned with (possibly NA)
  terms; if `"none"`, terms are not added.

- legacy:

  Logical scalar indicating whether to pull data from ExperimentHub. By
  default, we use data from the gypsum backend.

- ref:

  Reference data for annotation. One of:

  - `"HPCA"`: Uses the *Human Primary Cell Atlas* from the `celldex`
    package.

  - A `SingleCellExperiment` object containing pre-labeled reference
    cells.

  Note: The placeholder option `"custom"` is not valid—users must
  provide an actual reference object.

  A numeric matrix of (usually normalized and log-transformed)
  expression values from a reference dataset, or a SummarizedExperiment
  object containing such a matrix; see
  [`SingleR::trainSingleR`](https://rdrr.io/pkg/SingleR/man/trainSingleR.html)
  for details. Alternatively, a list or List of SummarizedExperiment
  objects or numeric matrices containing multiple references. Row names
  may be different across entries but only the intersection will be
  used,

- labels:

  A character vector of cell type labels corresponding to the columns of
  `ref`. Required when `ref` is a custom `SingleCellExperiment`. For
  `"HPCA"`, defaults to `ref$label.main`.

  Alternatively, if ref is a list, labels should be a list of the same
  length. Each element should contain a character vector or factor
  specifying the labels for the columns of the corresponding element of
  ref.

- ...:

  Additional arguments passed to
  [`SingleR::SingleR()`](https://rdrr.io/pkg/SingleR/man/SingleR.html)
  (e.g., `assay.type.ref`, `genes`).

## Value

- If input is a `Seurat` object: returns the same object with three new
  metadata columns:

  - `SingleR_labels`:

    Predicted cell type labels.

  - `SingleR_delta_next`:

    Confidence score (difference between top and second-best scores).

  - `SingleR_pruned_labels`:

    Labels after pruning low-confidence assignments.

- If input is a `SingleCellExperiment`: returns the `SingleR` prediction
  result

## Requirements

- The `celldex` package is required when using `ref = "HPCA"`.

- All gene names must be consistent between test and reference datasets
  (e.g., Ensembl or symbol).

## See also

Other Single_Cell_Annotation_Method:
[`CellTypistAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/CellTypistAnnotate.md),
[`RegisterAnnoMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterAnnoMethod.md),
[`SCAnnotateStrategy`](https://wanglabcsu.github.io/sigbridger/reference/SCAnnotateStrategy.md),
[`mLLMCelltypeAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/mLLMCelltypeAnnotate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Using HPCA reference (requires celldex)
annotated_seurat <- SingleRAnnotate(seurat_obj, ref = "HPCA")

# Using custom reference
library(celldex)
ref_sce <- BlueprintEncodeData()
annotated_sce <- SingleRAnnotate(
  sc = my_sce,
  ref = ref_sce,
  labels = ref_sce$label.main,
  de.n = 30
)
} # }
```
