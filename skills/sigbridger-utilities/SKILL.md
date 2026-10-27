---
name: sigbridger-utilities
description: Use SigBridgeR's helper toolkit around the screening core — LoadRefData for example data, SCIntegrate for merging heterogeneous datasets, ChooseNormalization and FindRobustElbow for preprocessing decisions, AddMisc/AddMetaFeature for Seurat metadata, bulk-expression helpers (SymbolConvert, AggregateDups, CountsToTPM, FPKMToTPM, GetGeneLength, CheckNA, IsCountsMatrix, DrawBulkPCA), PhenoMap for phenotype recoding, SCAnnotate for cell-type annotation, and global options via getFuncOption/setFuncOption/setThreads. Use when a SigBridgeR task needs data preparation, ID conversion, or runtime configuration rather than screening itself.
---

# SigBridgeR Utilities

## Purpose

These functions prepare inputs and configure the runtime. They are independent of
`Screen()` and can be used on their own.

## Reference data — `LoadRefData()`

```r
LoadRefData(data_type = c("survival", "binary", "continuous"),
            path = NULL, timeout = <option>)
```

Returns a list of three: `[[1]]` Seurat object, `[[2]]` bulk TPM matrix,
`[[3]]` phenotype (named vector, or `data.frame(time, status)` for survival).
Downloads once into a user cache (`~/.cache/SigBridgeR`, or
`%LOCALAPPDATA%\SigBridgeR\cache`), with GitHub Raw → GitHub API → jsDelivr
fallback. Use it to check the expected data shapes.

## Runtime options

```r
getFuncOption("verbose")        # read a package option
setFuncOption(verbose = FALSE)  # set package options globally
setThreads(threads = 8L)        # data.table / cheapr / qs2 / OpenMP / WGCNA / TF
```

`setThreads()` applies one global count to every backend unless overridden per
backend (`dt`, `cheapr`, `qs2`, `openmp`, `WGCNA`, `tf_config`). Call it **before**
importing TensorFlow through reticulate.

## Single-cell preparation

### Integrating heterogeneous datasets — `SCIntegrate()`

```r
integrated <- SCIntegrate(seu1, seu2, seu3)                      # Harmony by default
integrated <- SCIntegrate(seu1, seu2, method = Seurat::CCAIntegration)
merged_mat <- SCIntegrate(mat1, mat2)                            # matrix method
```

- Seurat inputs are merged, then processed by the `pipeline` letter string
  (default `"nsvpiectu"` = Normalize → Scale → VariableFeatures → PCA →
  IntegrateLayers → FindNeighbors → FindClusters → tSNE → UMAP).
- `method` defaults to `Seurat::HarmonyIntegration`; `add.cell.ids` prefixes are
  inferred from argument names.
- Matrix/data.frame inputs are gene-wise unioned and column-concatenated (missing
  values `NA`).
- Use this — not `MergeResult()` — when samples/cells differ between datasets.

### Normalisation and dimensionality

```r
res <- ChooseNormalization(sct = sct, lognorm = lognorm,
                           n_hvgs = 2000L, low_expressed_thresh = 0.2,
                           weight = c(variance_stability = 0.4,
                                      marker_signal = 0.35,
                                      dropout_robustness = 0.25))
res$metrics; res$recommendation

ndims <- FindRobustElbow(seurat_obj, ndims = 50L)   # requires RunPCA() first
```

`ChooseNormalization()` scores any number of equivalently preprocessed Seurat
objects (arguments are named per method) on variance stabilisation, marker
retention and dropout robustness; `weight` must sum to 1.

### Seurat bookkeeping

```r
seurat_obj <- AddMisc(seurat_obj, markers = c("CD3D", "CD4"))          # @misc
seurat_obj <- AddMisc(seurat_obj, markers = c("CD3D"), cover = FALSE)  # auto-suffix
seurat_obj <- AddMetaFeature(seurat_obj, gene_type = gene_type)        # assay-level
```

`AddMisc()` writes to `@misc` (supports multiple attributes at once, optional
auto-increment when `cover = FALSE`); `AddMetaFeature()` writes gene-level
metadata to `@assays[[assay]]@meta.data`, accepting vectors or 2-D objects
aligned by feature name. SigBridgeR itself stores screening evidence in `@misc`.

### QC helpers

```r
obj <- QCPatternDetect(obj, pattern = c("^MT-", "^RP[SL]"))
obj <- QCFilter(obj, data_filter.thresh = list(percent.mt = 20L, percent.rp = 60L))
Pattern2Colname("^MT-|^RP[LS]")     # "mt_rp" — how column names are generated
```

`SCPreProcess()` (the old pipeline engine) is **deprecated since 4.0.0** and emits
a lifecycle warning; use Seurat + `scCustomize` instead.

## Bulk expression helpers

| Function | Purpose |
| --- | --- |
| `SymbolConvert(data, genome_build = c("hg38","hg19","mm10","mm9","ce11","T2T"))` | Ensembl/TCGA IDs → gene symbols; unmatched IDs are kept, duplicates collapsed. Needs `IDConverter`. |
| `AggregateDups(x, method = c("max","sum","mean","median","first"))` | collapse duplicate rows/columns (`row_method`, `col_method` separately) |
| `CheckNA(data, max_print = 5L)` | report missing values |
| `IsCountsMatrix(x, integer_tol = 1e-8, min_integer_fraction = 0.95)` | detect raw counts vs. normalised values |
| `GetGeneLength(x, canonical_chr_only = FALSE)` | gene lengths for TPM conversion |
| `CountsToTPM(counts, gene_length)` | counts → TPM |
| `FPKMToTPM(data, na_as_zero = TRUE)` | FPKM → TPM |
| `DrawBulkPCA(bulk, group, batch = NULL)` | PCA scatter of bulk samples |

Run `SymbolConvert()` before `Screen()` whenever bulk and single-cell identifiers
differ — the screening step assumes shared gene symbols.

## Phenotype helpers — `PhenoMap()`

`dplyr::case_when`-style recoding, vectorised on data.tables:

```r
pheno <- PhenoMap(pheno, value > 60 ~ "High", value > 40 ~ "Medium",
                  .default = "Low")
df <- PhenoMap(df, age < 30 ~ "Young", age < 50 ~ "Middle", .default = "Senior")
```

Conditions must be formulas; the first match wins.

## Cell type annotation — `SCAnnotate()`

```r
SCAnnotate(sc, method = c("CellTypist", "SingleR", "mLLMCelltype"), ...)
names(SCAnnotateStrategy)     # registered methods
```

Dispatches to `SingleRAnnotate()`, `CellTypistAnnotate()` or
`mLLMCelltypeAnnotate()`. Result columns: `SingleR_labels` /
`SingleR_delta_next` / `SingleR_pruned_labels`; `celltypist_prediction` (+
probability); `mllmcelltype_cell_type` / `_consensus_proportion` / `_entropy`.
Partial method names are accepted. Python-backed methods need an environment:
`SetupPyEnv()`, `ListPyEnv()`, `DownloadCellTypistModels()`,
`ListCellTypistModels()`. Add your own method with `Register()` — see
`sigbridger-extend`.

## Order of operations for a real analysis

1. `SymbolConvert()` on bulk; check `IsCountsMatrix()`.
2. `QCPatternDetect()` + `QCFilter()`, `ChooseNormalization()`,
   `FindRobustElbow()`, then cluster/UMAP.
3. `PhenoMap()` to shape the phenotype; verify IDs against `colnames(bulk)`.
4. `setThreads()`, then `Screen()` (see `sigbridger-quickstart`).
5. `SCAnnotate()` if cell types are needed to interpret the screened fractions.

## Reference

- `?SigBridgeR::LoadRefData`, `?SigBridgeR::SCIntegrate`,
  `?SigBridgeR::ChooseNormalization`, `?SigBridgeR::FindRobustElbow`,
  `?SigBridgeR::AddMisc`, `?SigBridgeR::AddMetaFeature`,
  `?SigBridgeR::SymbolConvert`, `?SigBridgeR::PhenoMap`,
  `?SigBridgeR::SCAnnotate`, `?SigBridgeR::setThreads`
- Vignette `Other_Function_Details` and `Full_Tutorial` §1
