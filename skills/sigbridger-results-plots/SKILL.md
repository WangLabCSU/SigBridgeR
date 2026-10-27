---
name: sigbridger-results-plots
description: Inspect and combine SigBridgeR screening output — read a ScreenMethodResult S7 object, merge results from several methods with MergeResult(), aggregate categorical labels with WeightedVote(), and plot them with ScreenFractionPlot(), ScreenUpset() and autoplot(). Use when the user has screening results and asks how to extract the Positive/Negative/Neutral labels, compare methods, or draw fraction/UpSet/UMAP figures.
---

# SigBridgeR Results & Visualisation

## What `Screen()` hands back

`Screen()` returns a `ScreenMethodResult` (S7, inherits `SigBridgeRBase`):

| Access | Content |
| --- | --- |
| `res@scRNA_data` | `Seurat` object; labels in `meta.data` |
| `res$extra` | list of method-specific extras (scores, model, stats) |
| `res@extra` | same, as an S7 property |
| `print(res)` | property dump (cli-formatted) |
| `summary(res)` / `format(res)` | compact overview |
| `as.list(res)` | named list of properties |
| `names(res)` / `length(res)` | property names / count |
| `res1 + res2` | shorthand for `MergeResult(res1@scRNA_data, res2@scRNA_data)` |

Labels live in a metadata column whose name equals `label_type` (default: the
method name, e.g. `scissor`, `scPAS`). Values are `"Positive"`, `"Negative"`,
`"Neutral"` (occasionally `"Other"`). Never hard-code the column name — list it:

```r
colnames(res@scRNA_data[[]])
table(res@scRNA_data$scissor)
```

## Merging several methods — `MergeResult()`

```r
merged <- MergeResult(res1, res2, res3)                       # plain merge
merged <- MergeResult(res1, res2, weights = c(Scissor = 3, scPP = 1),
                      ties.method = "first")                  # with voting
```

Accepted inputs: `ScreenMethodResult`, `Seurat`, or a list containing
`scRNA_data`. Invalid inputs are skipped with a warning.

Behaviour to remember:

- **Inner join on cell barcodes** — only cells present in *all* inputs survive.
  Retain-count arithmetic: `ncol(merged) <= ncol(first_object)`.
- Expression data, assays, reductions, graphs and images come from the **first**
  object; unique slots from the others are added.
- Duplicate metadata columns keep the **first occurrence**.
- A warning about "not designed for integrating heterogeneous datasets" means the
  inputs have different cell sets — use `SCIntegrate()` (see
  `sigbridger-utilities`) for different samples, not `MergeResult()`.
- With `weights`, a `vote` column is appended with the aggregated label.

## Voting — `WeightedVote()`

```r
vote_data <- data.frame(
  scissor = c("Positive", "Negative"),
  scPP    = c("Positive", "Positive"),
  row.names = c("cell1", "cell2")
)
WeightedVote(vote_data, weights = c(scissor = 3, scPP = 1), ties.method = "first")
# -> "Positive" "Positive" (names = cell ids)
```

- Every cell value must be one of `"Positive"`, `"Negative"`, `"Neutral"`,
  `"Other"`.
- `names(weights)` must be a permutation of `colnames(vote_data)`; otherwise it
  aborts.
- `ties.method` ∈ `"random"`, `"first"`, `"last"`; the winner is reported as
  `"Positive"` or `"Other"`.
- Row names are preserved as names on the returned character vector.

`MergeResult(..., weights = )` calls this internally after matching metadata
column names against `ScreenStrategy` keys (case-insensitive).

## Plots

### Stacked fractions — `ScreenFractionPlot()`

```r
fraction <- ScreenFractionPlot(
  merged,
  group_by    = "seurat_clusters",     # any metadata column; default "Source"
  screen_type = c("scissor", "scPAS"), # exact metadata column names
  ncol        = 2,
  show_plot   = FALSE
)
names(fraction)       # stats, plot, combined_plot
fraction$combined_plot
```

Returns `stats` (cell counts and percentages per group), `plot` (one ggplot or a
list) and `combined_plot` (patchwork, multiple screens only). Bars are ordered by
descending Positive fraction; zero-fraction groups are dropped unless
`show_null = TRUE`. Requires `ggplot2`, `patchwork`, `dplyr`.

### UpSet intersections — `ScreenUpset()`

```r
upset <- ScreenUpset(merged, screen_type = c("scissor", "scPAS"))
names(upset)   # plot, stats
upset$plot
```

Computes every combination of the requested screens and counts cells positive for
all of them. `screen_type = NULL` auto-detects columns matching
`names(ScreenStrategy)`. Requires `ggupset`, `tibble`, `ggplot2`. Tune with
`n_intersections`, `order_by`, `bar_color`, `combmatrix_point_color` plus ggplot2
theme arguments.

### UMAP — `autoplot()`

```r
ggplot2::autoplot(merged, group.by = "scissor")   # Seurat method
ggplot2::autoplot(res)                            # ScreenMethodResult method
```

Wraps `Seurat::DimPlot()` with `tidydr::theme_dr()`. When `group.by` matches a
registered screening method the palette is forced to grey / `#c24b4b` (Positive) /
`#5189bb` (Negative); otherwise `palette_SigBridgeR()` is used. Requires `tidydr`.

Venn diagrams are not built in — extract positive cells and use
`ggVennDiagram`:

```r
pos <- lapply(c("scissor", "scPAS"),
              \(m) colnames(merged)[merged[[m]] == "Positive"])
ggVennDiagram::ggVennDiagram(pos)
```

## Export

`merged` is a plain `Seurat` object, so anything Seurat/scCustomize offers works:

```r
SeuratObject::SaveSeuratRds(merged, "screened_seurat.rds")
```

## Pitfalls

- Plotting the `ScreenMethodResult` directly with Seurat functions fails — pass
  `res@scRNA_data` or use `autoplot()`.
- `screen_type` values must be **exact** metadata column names (case included).
- `MergeResult()` silently reduces the cell set; assert
  `ncol(merged) == ncol(base)` when all methods ran on the same data.
- `weights` for methods that did not run are ignored with a warning.

## Reference

- `?SigBridgeR::MergeResult`, `?SigBridgeR::WeightedVote`,
  `?SigBridgeR::ScreenFractionPlot`, `?SigBridgeR::ScreenUpset`,
  `?SigBridgeR::autoplot`
- Vignette `Quick_Start` (screening → merge → visualise) and `Full_Tutorial`
