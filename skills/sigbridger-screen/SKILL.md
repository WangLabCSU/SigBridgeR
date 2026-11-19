---
name: sigbridger-screen
description: Run a SigBridgeR phenotype-associated cell screening analysis after validating its three required inputs. Check that single-cell data is a Seurat object, bulk RNA-seq is a genes-by-samples matrix, and phenotype is a named 0/1 vector (binary) or a data.frame whose first two columns are time and status (survival); verify the chosen method is in names(SigBridgeR::ScreenStrategy), call Screen(), and read the labelled Seurat object from the result's @scRNA_data slot. Use when the user wants to run Screen(), asks whether their scRNA/bulk/phenotype data are ready for screening, or asks which screening algorithms are available.
---

# SigBridgeR Screening

## Purpose

Screen phenotype-associated cells from single-cell RNA-seq by bridging a matched
bulk RNA-seq matrix and sample-level phenotype. The contract is strict: wrong data
shapes fail before any computation. This skill validates all three inputs, selects
a registered algorithm, runs `Screen()`, and locates the output in the
`ScreenMethodResult` object.

## When to use

- "Run screening / `Screen()` on my scRNA + bulk + phenotype data."
- "Are my data in the right shape for SigBridgeR?"
- "Which screening methods can I use?"

## Step 0 — Confirm all three inputs exist

Before anything else, confirm the session contains the single-cell object, the
matched bulk matrix and the phenotype vector/table. If one is absent, stop and ask
the user to provide it; do not guess or fabricate data.

```r
needed <- c("sc_data", "matched_bulk", "phenotype")
absent <- needed[!vapply(needed, exists, logical(1L), envir = .GlobalEnv)]
if (length(absent)) {
  rlang::abort("Missing screening input(s): {paste(absent, collapse = ', ')}")
}
```

## Step 1 — Validate the shapes

### 1a. Single-cell data must be a Seurat object

```r
if (!inherits(sc_data, "Seurat")) {
  rlang::abort("`sc_data` must be a Seurat object (class `Seurat`).")
}
```

Use a preprocessed object (QC, normalisation, clustering, UMAP already done);
`Screen()` does not build those for you.

### 1b. Bulk RNA-seq must be a genes x samples matrix

Rows are genes, columns are samples. Both dimensions need names: row names are
matched against single-cell features, column names against phenotype names.

```r
if (!is.matrix(matched_bulk)) {
  rlang::abort("`matched_bulk` must be a matrix with genes in rows, samples in columns.")
}
if (is.null(rownames(matched_bulk)) || is.null(colnames(matched_bulk)) ||
    nrow(matched_bulk) == 0L || ncol(matched_bulk) == 0L) {
  rlang::abort("`matched_bulk` needs non-empty gene row names and sample column names.")
}
```

### 1c. Phenotype depends on `phenotype_class`

**Binary (`phenotype_class = "binary"`) — named vector encoded as 0/1:**

```r
if (!is.numeric(phenotype) || !is.vector(phenotype) || is.null(names(phenotype))) {
  rlang::abort("Binary `phenotype` must be a named numeric vector with values 0/1.")
}
if (!setequal(names(phenotype), colnames(matched_bulk))) {
  rlang::abort("`names(phenotype)` must match `colnames(matched_bulk)` (same samples).")
}
if (!all(phenotype %in% c(0L, 1L)) || length(unique(phenotype)) != 2L) {
  rlang::abort("Binary `phenotype` must contain only 0 and 1, with both groups present.")
}
phenotype <- phenotype[colnames(matched_bulk)]   # align to bulk column order
```

**Survival (`phenotype_class = "survival"`) — data.frame with `time`, `status` first:**

```r
if (!is.data.frame(phenotype) || ncol(phenotype) < 2L ||
    !identical(names(phenotype)[1:2], c("time", "status"))) {
  rlang::abort("Survival `phenotype` must be a data.frame whose first two columns are `time` and `status`.")
}
if (!setequal(rownames(phenotype), colnames(matched_bulk))) {
  rlang::abort("`rownames(phenotype)` must match `colnames(matched_bulk)` (same samples).")
}
phenotype <- phenotype[colnames(matched_bulk), , drop = FALSE]
```

`status` uses 0 = censored, 1 = event. Other phenotype classes (`continuous`)
use a named numeric vector, but are not the scope of this skill.

### 1d. Check genes and samples are genuinely matched

```r
features <- rownames(Seurat::GetAssayData(sc_data, layer = "data"))
if (length(intersect(rownames(matched_bulk), features)) == 0L) {
  rlang::abort("No shared gene IDs between bulk and single-cell data; run `SymbolConvert()` first.")
}
```

## Step 2 — Select an algorithm

`screen_method` must be a key currently registered in `ScreenStrategy`. Query it at
runtime rather than hard-coding a name — `Register()` may add methods:

```r
methods <- names(SigBridgeR::ScreenStrategy)
if (!screen_method %in% methods) {
  rlang::abort("`screen_method` must be one of: {paste(methods, collapse = ', ')}")
}
cfg <- SigBridgeR::ScreenStrategy[[screen_method]]
if (!phenotype_class %in% cfg@phenotype_class) {
  rlang::abort("{screen_method} does not support phenotype_class = '{phenotype_class}'.")
}
```

Built-ins at the time of writing: `Scissor`, `scPAS`, `scAB`, `scPP`, `DEGAS`,
`DEGASv2`, `LP_SGL`, `PIPET`, `SIDISH`, `SCIPAC`, `TiRank`. Restrictions exist
(`scAB`: binary/survival only; `SIDISH`: survival only), so always use the check
above. Adding a new method is a different workflow — use the `sigbridger-extend`
skill.

## Step 3 — Run `Screen()`

```r
res <- SigBridgeR::Screen(
  matched_bulk    = matched_bulk,
  sc_data         = sc_data,
  phenotype       = phenotype,
  phenotype_class = phenotype_class,   # "binary" or "survival"
  screen_method   = screen_method,
  ...                                   # method-specific parameters
)
```

Always use named arguments. `label_type` defaults to `screen_method`; set it
explicitly when running the same method on different phenotypes, otherwise the
metadata label column is overwritten.

## Step 4 — Read the result from `@scRNA_data`

`Screen()` returns a `ScreenMethodResult` S7 object, **not** a bare Seurat object.
The screened, labelled Seurat object lives in the `scRNA_data` slot; method-specific
extras live in `$extra`.

```r
screened <- res@scRNA_data
labels <- if (is.null(label_type)) screen_method else label_type
colnames(screened[[]])          # find the exact label column
table(screened[[labels]])       # Positive / Negative / Neutral
res$extra                       # scores, model, stats (method-dependent)
```

Pass `res@scRNA_data` (or use `autoplot(res)`) to Seurat plotting functions —
calling them on `res` directly fails.

## Failure signals

| Signal | Cause | Fix |
| --- | --- | --- |
| `class` check fails | object is not Seurat | preprocess with Seurat first |
| `value must be a matrix` | bulk is a data.frame | `as.matrix()` with gene/sample names |
| names mismatch | sample IDs differ | reorder by `colnames(matched_bulk)` as above |
| `phenotype_class not supported` | method/trait mismatch | choose another `screen_method` |
| no shared gene IDs | Ensembl vs symbol | `SymbolConvert()` on bulk or single-cell features |

## Reference

- Source: `R/20-QuickScreen.R` (`Screen()`), `R/70-ScreenStrategy.R` (`ScreenStrategy`),
  `R/04e-ScreenMethodResult.R` (`scRNA_data` slot), `R/04a-ClassBase.R` (input validators).
- In-session: `?SigBridgeR::Screen`, `?SigBridgeR::ScreenStrategy`,
  `?SigBridgeR::ScreenMethodResult`, `?SigBridgeR::SymbolConvert`.
- Vignettes: `Quick_Start`, `Full_Tutorial`.
