    library(SigBridgeR)
    library(Seurat)

# Overview

`BenchGeneScore()` scores one or more gene sets in **every cell** of a
screened single-cell dataset. It does not compare labels against a
ground truth like `BenchAcc()`; instead it asks whether the cells
selected by different methods carry the pathway activity expected from
the phenotype.

Four scoring backends are available:

<table>
<colgroup>
<col style="width: 38%" />
<col style="width: 34%" />
<col style="width: 26%" />
</colgroup>
<thead>
<tr>
<th><code>method</code></th>
<th>Backend</th>
<th>Notes</th>
</tr>
</thead>
<tbody>
<tr>
<td><code>"ssGSEA"</code></td>
<td><code>GSVA::ssgseaParam()</code> + <code>GSVA::gsva()</code></td>
<td>Rank-based, returns gene sets by cells</td>
</tr>
<tr>
<td><code>"AUCell"</code></td>
<td><code>AUCell::AUCell_buildRankings()</code> +
<code>AUCell::AUCell_calcAUC()</code></td>
<td>AUC of the gene-set ranking</td>
</tr>
<tr>
<td><code>"UCell"</code></td>
<td><code>UCell::ScoreSignatures_UCell()</code></td>
<td>Fast rank-based scores, transposed to the common layout</td>
</tr>
<tr>
<td><code>"AddModuleScore"</code></td>
<td><code>Seurat::AddModuleScore()</code></td>
<td>Control-gene corrected mean expression</td>
</tr>
</tbody>
</table>

All backends return the same orientation: gene sets in rows and cells in
columns. That makes it possible to switch backend or compare methods
without reshaping the result.

# A small example

The following object has twenty genes and eight cells. The first pathway
is up-regulated in the `"Positive"` cells, while the second pathway is
not, so the score matrix should separate the two groups for the first
pathway.

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

`BenchGeneScore()` accepts a `Seurat` object, a `ScreenMethodResult`, or
a `BenchGeneScoreConfig`. The `pathway` argument is a named list of gene
sets:

    scores <- BenchGeneScore(
      sc,
      pathway = pathways,
      method = "AddModuleScore",
      nbin = 4L,
      ctrl = 2L,
      verbose = FALSE
    )

    round(scores, 3)

The result is a gene-set-by-cell matrix whose columns match the cells of
`sc`:

    identical(colnames(scores), colnames(sc))

The group means show the expected pattern: `Pathway_A` is higher in
`"Positive"` cells, while `Pathway_B` does not separate the groups.

    round(
      apply(scores, 1L, function(x) tapply(x, sc$truth, mean)),
      3
    )

# Choosing a backend

The backend is selected with `method`. Additional named arguments in
`...` are matched against the formal arguments of the backend and, for
two-stage backends, split over both stages automatically:

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

The required package is loaded on demand with
`rlang::check_installed()`, so the generic can be called even when only
one backend is installed. In an interactive session, choose the backend
that matches the question and keep it fixed when comparing screening
methods:

- `"ssGSEA"` and `"UCell"` are rank-based and robust to library-size and
  technical variation. They are a good default for comparing cell
  populations.
- `"AUCell"` also uses rankings and is useful when the gene-set signal
  is concentrated in the top-expressed genes.
- `"AddModuleScore"` is convenient when the object is already processed
  with Seurat, but it is sensitive to the choice of control genes and
  binning parameters (`nbin`, `ctrl`).

# Configuration input

A `BenchGeneScoreConfig` runs the configured screening methods first and
scores the accumulated labels in a single pass, so the scores are
comparable across methods by construction:

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

`BenchGeneScoreConfigSet()` stores the backend and its arguments; the
scoring backend’s package is checked when the configuration is built, so
a missing dependency is reported before the screening run starts.

# Reading the scores

A gene-set-by-cell matrix is deliberately simple: each value is a
relative activity score for one gene set in one cell. It is not a
p-value and should not be compared across datasets or across different
backends as if it were on a common absolute scale.

Useful summaries are:

    # * mean score per ground-truth or cluster label
    cluster_mean <- apply(
      scores,
      1L,
      function(x) tapply(x, sc$truth, mean)
    )
    round(cluster_mean, 3)

    # * rank cells within each gene set
    pathway_rank <- t(apply(scores, 1L, rank))
    pathway_rank[, 1:4]

A common workflow is to compare the distribution of a pathway score
across the cells selected by two methods. Because the matrix shares its
columns with the `Seurat` object, the labels can be attached directly:

    sc$Pathway_A <- scores["Pathway_A", ]

    aggregate(
      Pathway_A ~ truth,
      data = sc[[]],
      FUN = mean
    )

# Practical notes

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

# See also

- `BenchGeneScoreConfig()` and `BenchGeneScoreConfigSet()` for the
  configuration interface.
- `BenchGSEA()` for pathway enrichment of method marker genes.
- `BenchAcc()` for classification accuracy against a ground truth.
- `BenchMemory()` for run time and memory.
