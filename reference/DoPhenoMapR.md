# PhenoMapR: Bulk-to-Single-Cell Phenotype Mapping

Derives a gene-level reference z-score from matched bulk expression and
a phenotype (binary, continuous, or survival), then scores single cells
and assigns tail phenotype groups via PhenoMapR.

## Usage

``` r
DoPhenoMapR(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = "PhenoMapR",
  phenotype_class = c("binary", "survival", "continuous"),
  phenotype_column = NULL,
  survival_time = "time",
  survival_event = "status",
  normalize = TRUE,
  platform = c("auto", "rnaseq", "microarray"),
  probe_annotation = NULL,
  hugo_species = c("human", "mouse"),
  binary_positive_reference = c("second", "first"),
  assay = "RNA",
  layer = "data",
  z_score_cutoff = 2L,
  reference_sign = 1L,
  score_mode = c("weighted_sum", "activity_adjusted"),
  vars_to_regress = c("S.Score", "G2M.Score", "nCount_RNA"),
  cancer_type = NULL,
  pseudobulk = FALSE,
  group_by = NULL,
  sample_metadata = NULL,
  percentile = 0.05,
  permutation_n = NULL,
  permutation_seed = NULL,
  ...,
  seed = SigBridgeRUtils::getFuncOption("seed"),
  verbose = SigBridgeRUtils::getFuncOption("verbose")
)
```

## Arguments

- matched_bulk:

  Matrix or data frame of bulk expression (genes x samples). Column
  names must match `phenotype` names/row names.

- sc_data:

  A Seurat object containing the cells to be scored.

- phenotype:

  Phenotype data. A named vector for binary/continuous traits, or a data
  frame with time-to-event and event columns for survival traits;
  names/row names must match `colnames(matched_bulk)`.

- label_type:

  Character. `meta.data` column receiving the phenotype group. Default:
  `"PhenoMapR"`.

- phenotype_class:

  One of `"binary"`, `"survival"`, or `"continuous"`.

- phenotype_column:

  Character. Column name in `phenotype` for the outcome. For
  `phenotype_type = "survival"` this is ignored; use `survival_time` and
  `survival_event` instead.

- survival_time:

  Character. Column name in `phenotype` for time-to-event (e.g. overall
  survival time). Required when `phenotype_type = "survival"`.

- survival_event:

  Character. Column name in `phenotype` for event indicator (0/1 or
  FALSE/TRUE). Required when `phenotype_type = "survival"`.

- normalize:

  Logical. If `TRUE`, run normalization when expression looks like
  counts (default `TRUE`). Set `FALSE` to skip.

- platform:

  One of `"auto"`, `"rnaseq"`, or `"microarray"`. Microarray inputs
  follow the PRECOG preprocessing strategy: optional probe-to-gene
  mapping, log2 transform when needed, removal of zero-variance genes
  and samples/genes with heavy missingness, per-study quantile
  normalization, and per-gene unit variance scaling across samples.
  RNA-seq inputs use log2(CPM+1) when the matrix looks like raw counts.

- probe_annotation:

  Optional data.frame with probe IDs and gene symbols (required for
  probe-level microarray matrices unless row names are already gene
  symbols). See
  [`collapse_probes_to_genes`](https://brooksbenard.github.io/PhenoMapR/reference/collapse_probes_to_genes.html).

- hugo_species:

  Character. Species for HUGO symbol cleaning: `"human"` or `"mouse"`
  (default `"human"`).

- binary_positive_reference:

  For `phenotype_type` `"binary"` only: which level of the binary factor
  should correspond to the **positive** outcome in logistic regression
  (`y = 1`), so that genes with **positive** z-scores are those whose
  **higher** expression is associated with that level. Use `"first"`
  when the first level of `factor(..., levels = c(...))` is your
  phenotype of interest (e.g. `mutated` vs `wt`); use `"second"` for the
  legacy convention (second level coded as `y = 1`). Default `"second"`
  preserves behaviour from previous versions when levels were implicit
  (e.g. alphabetical).

- assay:

  Assay name for Seurat/SCE objects (default: "RNA" for sc, "Spatial"
  for spatial)

- layer:

  Seurat matrix to score (`"data"` normalized, `"counts"` raw, or
  `"scale.data"`; default `"data"`). Seurat v5 terminology. PhenoMapR
  maps this to `GetAssayData(layer=...)` or `GetAssayData(slot=...)`
  based on the installed SeuratObject version.

- z_score_cutoff:

  Absolute z-score threshold for filtering genes (default: 2). Built-in
  PRECOG/TCGA matrices shipped in the package are already sparse at this
  cutoff (entries with \\\|z\| \leq 2\\ stored as `NA`; see
  [`precog`](https://brooksbenard.github.io/PhenoMapR/reference/precog.html),
  [`tcga`](https://brooksbenard.github.io/PhenoMapR/reference/tcga.html)),
  so lowering the cutoff below 2 does not recover discarded associations
  unless a denser custom reference is supplied.

- reference_sign:

  Multiply all reference gene z-scores by this value before the weighted
  sum. Use `1` (default) or `-1` to flip the sign of resulting PhenoMap
  scores without changing the reference matrix (e.g. if you want
  positive scores to mean the opposite association). For custom
  references from
  [`derive_reference_from_bulk`](https://brooksbenard.github.io/PhenoMapR/reference/derive_reference_from_bulk.html),
  prefer setting `binary_positive_reference` so the reference itself
  matches your factor levels (first level = positive association when
  `"first"`).

- score_mode:

  Scoring mode: `"weighted_sum"` (default) uses the standard dot product
  of expression and reference z-scores; `"activity_adjusted"` regresses
  technical covariates per gene (via Seurat `ScaleData`) before
  computing the weighted sum on the regressed layer. Requires Seurat for
  matrix inputs.

- vars_to_regress:

  Metadata columns to regress when `score_mode = "activity_adjusted"`.
  Default: `c("S.Score", "G2M.Score", "nCount_RNA")`. For matrix input,
  pass covariates via `sample_metadata` (see below). `nCount_RNA` is
  computed from column sums when absent.

- cancer_type:

  Cancer type label. Required if using built-in reference datasets.
  Should match column names in reference data. For ICI PRECOG, use
  format "ABBREV" or "ABBREV_Metastatic" (e.g., "MELANOMA",
  "MELANOMA_Metastatic")

- pseudobulk:

  Logical. If TRUE, aggregate expression before scoring (default: FALSE)

- group_by:

  Column name for pseudobulk grouping. Required if pseudobulk = TRUE.
  For Seurat/SCE objects, should be a metadata column.

- sample_metadata:

  Optional data.frame of per-cell/sample covariates for
  `score_mode = "activity_adjusted"` when `expression` is a matrix. Must
  include a cell identifier column (default `"cell_id"`).

- percentile:

  Fraction of cells in each tail (default 0.05 for top and bottom 5%).

- permutation_n:

  If a positive integer, compute gene-shuffle permutation null scores
  and add `empirical_p_*` columns alongside raw score columns.

- permutation_seed:

  Random seed for permutation null (default `42`).

- ...:

  Additional arguments reserved for method compatibility.

- seed:

  Integer. Random seed. Default: `getFuncOption("seed")`.

- verbose:

  Logical. Whether to print progress messages. Default:
  `getFuncOption("verbose")`.

## Value

A list containing:

- `scRNA_data`: Seurat object with `label_type` (`"Positive"`,
  `"Negative"`, `"Neutral"`, or `NA`) and a `PhenoMapR_*` score column
  in `meta.data`.

- `phenomap_scores`: Per-cell scores from PhenoMap.

- `reference`: Gene-level reference z-scores.

## Details

`define_phenotype_groups()` labels the upper tail of the score
distribution `"Most Adverse"` and the lower tail `"Most Favorable"`. The
score increases towards the phenotype-valued-as-1 end when the reference
is derived with the default orientation (for a numeric `0`/`1`
phenotype, `binary_positive_reference = "second"` regresses on the `"1"`
level, and for survival traits the Cox z-score is positive for
hazard-increasing genes). Therefore, when phenotype `1` denotes the
unfavourable direction, the `"Most Adverse"` tail is the one to be
labelled `"Positive"`; the group labels are mapped to the package-wide
vocabulary as `"Most Adverse" -> "Positive"`,
`"Most Favorable" -> "Negative"` and `"Other" -> "Neutral"`. Setting
`reference_sign = -1` mirrors the score axis, in which case the two
tails swap. If a phenotype is encoded so that the adverse class is *not*
the reference-positive level (e.g. `binary_positive_reference = "first"`
with numeric `0`/`1` labels), the score direction is reversed as well;
choose the arguments so that higher scores track the adverse phenotype.

## References

PhenoMapR: mapping phenotypic information from bulk to single-cell data.
<https://github.com/brooksbenard/PhenoMapR>

## See also

Other screen_method:
[`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md),
[`DoLP_SGL()`](https://wanglabcsu.github.io/sigbridger/reference/DoLP_SGL.md),
[`DoPIPET()`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md),
[`DoSCIPAC()`](https://wanglabcsu.github.io/sigbridger/reference/DoSCIPAC.md),
[`DoSIDISH()`](https://wanglabcsu.github.io/sigbridger/reference/DoSIDISH.md),
[`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md),
[`DoTiRank()`](https://wanglabcsu.github.io/sigbridger/reference/DoTiRank.md),
[`DoscAB()`](https://wanglabcsu.github.io/sigbridger/reference/DoscAB.md),
[`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md),
[`DoscPP()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPP.md)

Other PhenoMapR:
[`ValidatePhenoMapRParams()`](https://wanglabcsu.github.io/sigbridger/reference/ValidatePhenoMapRParams.md)

## Examples

``` r
if (FALSE) { # \dontrun{
res <- DoPhenoMapR(
  matched_bulk = bulk_matrix,
  sc_data = seurat_obj,
  phenotype = a_named_vector,
  label_type = "PhenoMapR",
  phenotype_class = "binary"
)
} # }
```
