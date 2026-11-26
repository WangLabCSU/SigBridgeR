# Benchmarking Screening Methods: Gene Set Scoring

``` r

library(SigBridgeR)
library(Seurat)
```

## Overview

[`BenchGeneScore()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScore.md)
scores one or more gene sets in **every cell** of a screened single-cell
dataset. It does not compare labels against a ground truth like
[`BenchAcc()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAcc.md);
instead it asks whether the cells selected by different methods carry
the pathway activity expected from the phenotype.

Four scoring backends are available:

| `method` | Backend | Notes |
|----|----|----|
| `"ssGSEA"` | [`GSVA::ssgseaParam()`](https://rdrr.io/pkg/GSVA/man/ssgseaParam-class.html) + [`GSVA::gsva()`](https://rdrr.io/pkg/GSVA/man/gsva.html) | Rank-based, returns gene sets by cells |
| `"AUCell"` | [`AUCell::AUCell_buildRankings()`](https://rdrr.io/pkg/AUCell/man/AUCell_buildRankings.html) + [`AUCell::AUCell_calcAUC()`](https://rdrr.io/pkg/AUCell/man/AUCell_calcAUC.html) | AUC of the gene-set ranking |
| `"UCell"` | [`UCell::ScoreSignatures_UCell()`](https://rdrr.io/pkg/UCell/man/ScoreSignatures_UCell.html) | Fast rank-based scores, transposed to the common layout |
| `"AddModuleScore"` | [`Seurat::AddModuleScore()`](https://satijalab.org/seurat/reference/AddModuleScore.html) | Control-gene corrected mean expression |

All backends return the same orientation: gene sets in rows and cells in
columns. That makes it possible to switch backend or compare methods
without reshaping the result.

## A small example

The following object has twenty genes and eight cells. The first pathway
is up-regulated in the `"Positive"` cells, while the second pathway is
not, so the score matrix should separate the two groups for the first
pathway.

``` r

truth <- rep(c("Positive", "Negative"), times = 4L)

counts <- Matrix::Matrix(
  matrix(
    rep(seq_len(20L), times = 8L) + rep(seq_len(8L), each = 20L),
    nrow = 20L,
    dimnames = list(paste0("Gene", seq_len(20L)), paste0("Cell", seq_len(8L)))
  ),
  sparse = TRUE
)

# * pathway A is active in the Positive cells
counts[1:5, truth == "Positive"] <- counts[1:5, truth == "Positive"] + 30L

sc <- CreateSeuratObject(counts = counts)
sc <- NormalizeData(sc, verbose = FALSE)
sc$truth <- truth

pathways <- list(
  Pathway_A = paste0("Gene", 1:5),
  Pathway_B = paste0("Gene", 6:10)
)
```

[`BenchGeneScore()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScore.md)
accepts a `Seurat` object, a `ScreenMethodResult`, or a
`BenchGeneScoreConfig`. The `pathway` argument is a named list of gene
sets:

``` r

scores <- BenchGeneScore(
  sc,
  pathway = pathways,
  method = "AddModuleScore",
  nbin = 4L,
  ctrl = 2L,
  verbose = FALSE
)

round(scores, 3)
#>            Cell1  Cell2  Cell3  Cell4  Cell5  Cell6  Cell7  Cell8
#> Pathway_A  0.734 -0.667  0.664 -0.533  0.607 -0.446  0.560 -0.384
#> Pathway_B -0.238  0.100 -0.211  0.077 -0.190  0.062 -0.173  0.053
```

The result is a gene-set-by-cell matrix whose columns match the cells of
`sc`:

``` r

identical(colnames(scores), colnames(sc))
#> [1] TRUE
```

The group means show the expected pattern: `Pathway_A` is higher in
`"Positive"` cells, while `Pathway_B` does not separate the groups.

``` r

round(
  apply(scores, 1L, function(x) tapply(x, sc$truth, mean)),
  3
)
#>          Pathway_A Pathway_B
#> Negative    -0.507     0.073
#> Positive     0.641    -0.203
```

## Choosing a backend

The backend is selected with `method`. Additional named arguments in
`...` are matched against the formal arguments of the backend and, for
two-stage backends, split over both stages automatically:

``` r

# * ssGSEA: extra arguments go to ssgseaParam() or gsva()
ssgsea_scores <- BenchGeneScore(
  sc,
  pathway = pathways,
  method = "ssGSEA",
  minSize = 5L,
  verbose = FALSE
)

# * AUCell: arguments are split between buildRankings() and calcAUC()
aucell_scores <- BenchGeneScore(
  sc,
  pathway = pathways,
  method = "AUCell",
  nCores = 1L,
  verbose = FALSE
)

# * UCell: arguments go to ScoreSignatures_UCell()
ucell_scores <- BenchGeneScore(
  sc,
  pathway = pathways,
  method = "UCell",
  maxRank = 500L,
  verbose = FALSE
)
```

The required package is loaded on demand with
[`rlang::check_installed()`](https://rlang.r-lib.org/reference/is_installed.html),
so the generic can be called even when only one backend is installed. In
an interactive session, choose the backend that matches the question and
keep it fixed when comparing screening methods:

- `"ssGSEA"` and `"UCell"` are rank-based and robust to library-size and
  technical variation. They are a good default for comparing cell
  populations.
- `"AUCell"` also uses rankings and is useful when the gene-set signal
  is concentrated in the top-expressed genes.
- `"AddModuleScore"` is convenient when the object is already processed
  with Seurat, but it is sensitive to the choice of control genes and
  binning parameters (`nbin`, `ctrl`).

## Configuration input

A `BenchGeneScoreConfig` runs the configured screening methods first and
scores the accumulated labels in a single pass, so the scores are
comparable across methods by construction:

``` r

config <- BenchGeneScoreConfig(
  matched_bulk = bulk,
  sc_data = sc,
  phenotype = phenotype,
  label_type = "{method_name}",
  phenotype_class = "binary",
  screen_method = c("Scissor", "scPAS"),
  method_param = MethodParamSet(
    Scissor_param = list(alpha = 0.05),
    scPAS_param = list(alpha = 0.01)
  ),
  bench_config = BenchGeneScoreConfigSet(
    method = "ssGSEA",
    minSize = 5L
  )
)

scores <- BenchGeneScore(config, verbose = FALSE)
```

[`BenchGeneScoreConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfigSet.md)
stores the backend and its arguments; the scoring backend’s package is
checked when the configuration is built, so a missing dependency is
reported before the screening run starts.

## Reading the scores

A gene-set-by-cell matrix is deliberately simple: each value is a
relative activity score for one gene set in one cell. It is not a
p-value and should not be compared across datasets or across different
backends as if it were on a common absolute scale.

Useful summaries are:

``` r

# * mean score per ground-truth or cluster label
cluster_mean <- apply(
  scores,
  1L,
  function(x) tapply(x, sc$truth, mean)
)
round(cluster_mean, 3)
#>          Pathway_A Pathway_B
#> Negative    -0.507     0.073
#> Positive     0.641    -0.203

# * rank cells within each gene set
pathway_rank <- t(apply(scores, 1L, rank))
pathway_rank[, 1:4]
#>           Cell1 Cell2 Cell3 Cell4
#> Pathway_A     8     1     7     2
#> Pathway_B     1     8     2     7
```

A common workflow is to compare the distribution of a pathway score
across the cells selected by two methods. Because the matrix shares its
columns with the `Seurat` object, the labels can be attached directly:

``` r

sc$Pathway_A <- scores["Pathway_A", ]

aggregate(
  Pathway_A ~ truth,
  data = sc[[]],
  FUN = mean
)
#>      truth  Pathway_A
#> 1 Negative -0.5074308
#> 2 Positive  0.6413492
```

## Practical notes

- **Keep the backend fixed.** ssGSEA, AUCell, UCell and AddModuleScore
  use different scales. Comparing methods is meaningful only when every
  method is scored by the same backend with the same arguments.
- **Match the identifier space.** Pathway genes must match the row names
  of the assay. Convert symbols/Ensembl IDs or select the matching
  `assay` and `layer` first.
- **Missing genes are handled differently.** `AddModuleScore` drops
  genes and gene sets that are absent from the assay; other backends may
  return `NA` or silently score fewer genes. Inspect the result whenever
  a pathway disappears.
- **The `Seurat` metadata is left untouched.** `AddModuleScore` is
  called with an internal column prefix, so the score matrix is returned
  without adding columns to `sc`.
- **Scores are relative.** High activity means “high within this dataset
  and this backend”, not “biologically active above a universal
  threshold”. Use group comparisons and ranks rather than absolute
  cut-offs.

## See also

- [`BenchGeneScoreConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfig.md)
  and
  [`BenchGeneScoreConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfigSet.md)
  for the configuration interface.
- [`BenchGSEA()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEA.md)
  for pathway enrichment of method marker genes.
- [`BenchAcc()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAcc.md)
  for classification accuracy against a ground truth.
- [`BenchMemory()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemory.md)
  for run time and memory.
