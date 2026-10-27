---
name: sigbridger-quickstart
description: Run a minimal end-to-end SigBridgeR workflow in R — load the bundled example single-cell/bulk/phenotype data, call Screen() with a method such as Scissor, inspect the ScreenMethodResult S7 object, merge results from several methods and visualise the screened cells. Use when the user asks to start with SigBridgeR, wants a reproducible minimal example, or asks how to run Screen() for the first time.
---

# SigBridgeR Quick Start

## Purpose

Produce a working, end-to-end SigBridgeR analysis: phenotype-associated cells are
screened from single-cell RNA-seq data by bridging bulk expression + sample-level
phenotype. This skill covers the shortest correct path. Deeper topics live in the
sibling skills `sigbridger-screen-methods`, `sigbridger-cache-config`,
`sigbridger-results-plots`, `sigbridger-extend`, `sigbridger-utilities`.

## When to use

- "How do I start with SigBridgeR?" / "give me a minimal example of Screen()".
- Reproducing the Quick Start vignette or debugging a first run.
- Verifying an installation.

## Prerequisites

```r
pak::pak("WangLabCSU/SigBridgeR")           # development version (recommended)
# or: install.packages("SigBridgeR", repos = "https://wanglabcsu.r-universe.dev")
```

Loading the package prints a startup message with the SigBridgeR version — that
message is the success signal. `library(Seurat)` is required for the returned
object.

## Workflow

### 1. Fetch reference data

`LoadRefData()` downloads a small example bundle to a local cache
(`~/.cache/SigBridgeR` on Unix, `%LOCALAPPDATA%\SigBridgeR\cache` on Windows)
and returns a list of exactly three elements:

```r
library(SigBridgeR)

ref        <- LoadRefData(data_type = "binary")  # "binary" | "survival" | "continuous"
seurat_obj <- ref[[1]]  # preprocessed Seurat object (QC, normalization, clusters, UMAP)
bulk_bi    <- ref[[2]]  # bulk TPM matrix: genes x samples
pheno_bi   <- ref[[3]]  # named vector (binary / continuous)
```

For `data_type = "survival"`, `ref[[3]]` is a `data.frame` with columns
`time` and `status` (0 = censored, 1 = event). Second call is offline — the file
is reused. Pass `path = "<dir>"` to control the cache location.

### 2. Run one screening method

```r
scissor_res <- Screen(
  matched_bulk     = bulk_bi,
  sc_data          = seurat_obj,
  phenotype        = pheno_bi,
  phenotype_class  = "binary",
  screen_method    = "Scissor",
  alpha            = 0.05
)
```

- Argument order is positional: `matched_bulk`, `sc_data`, `phenotype`. Prefer
  named arguments — it prevents swapping bulk and phenotype.
- `phenotype_class` must be one the method supports (see
  `sigbridger-screen-methods`).
- `label_type` defaults to `screen_method`; it names the metadata column and the
  `@misc` entry, so set it explicitly when running the same method twice on
  different phenotypes.

### 3. Inspect the result

`Screen()` returns a `ScreenMethodResult` S7 object, not a bare Seurat object.

```r
class(scissor_res)          # "ScreenMethodResult" ...
print(scissor_res)          # property-by-property dump
scissor_res@scRNA_data      # the Seurat object carrying the labels
scissor_res$extra           # method-specific extras captured from the method
table(scissor_res@scRNA_data$scissor)   # Positive / Negative / Neutral
```

`$`, `[[`, `names()`, `length()`, `as.list()` all work on the result object.

### 4. Run a second method and merge

```r
scpas_res <- Screen(
  matched_bulk = bulk_bi, sc_data = seurat_obj, phenotype = pheno_bi,
  phenotype_class = "binary", screen_method = "scPAS"
)

merged_seurat <- MergeResult(scissor_res, scpas_res)
# equivalent shorthand for two results:
# merged_seurat <- scissor_res + scpas_res
```

`MergeResult()` performs an **inner join on cell barcodes** and keeps the
expression data of the first input. It returns a plain `Seurat` object.

### 5. Visualise

```r
fraction <- ScreenFractionPlot(
  merged_seurat,
  group_by    = "seurat_clusters",
  screen_type = c("scissor", "scPAS")
)
fraction$combined_plot

upset <- ScreenUpset(merged_seurat, screen_type = c("scissor", "scPAS"))
upset$plot

ggplot2::autoplot(merged_seurat, group.by = "scPAS")   # UMAP coloured by label
```

`screen_type` must equal the metadata column name exactly. Do not guess it: run
`colnames(merged_seurat[[]])` and pass the label column(s) as written (typically
the `label_type` / method name, e.g. `"scissor"`, `"scPAS"`; some methods emit
extra columns such as `<label>_Positive`).

## Minimal checklist

1. `matched_bulk` columns == names/row names of `phenotype` (mismatch aborts
   before any computation).
2. Single-cell data is already normalised and clustered.
3. Gene identifiers agree between bulk and single-cell data; convert them with
   `SymbolConvert()` first if they are Ensembl/TCGA IDs.
4. Long-running methods: set `setThreads()` before running, and use `save_cache`
   so a failed run can be resumed.

## Common failure signals

| Message | Cause | Fix |
| --- | --- | --- |
| mismatch between bulk columns and phenotype names | sample IDs differ | reorder/rename before `Screen()` |
| `phenotype_class` not supported | e.g. `continuous` with `scAB` | pick another method |
| long silent run | algorithm iteration | run with `verbose = TRUE`, or cache |
| `ScreenMethodResult` has no metadata column | `label_type` overwritten | read `names(result$extra)` / `@misc` |

## Reference

- Package site: <https://wanglabcsu.github.io/SigBridgeR/>
- Vignettes: `Quick_Start`, `Full_Tutorial`, `Optimal_Params`, `Troubleshooting`
- In-session help: `?SigBridgeR::Screen`, `?SigBridgeR::MergeResult`
