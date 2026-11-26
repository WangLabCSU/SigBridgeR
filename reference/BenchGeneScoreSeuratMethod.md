# Gene Set Scoring Backends

Internal implementations behind
[`BenchGeneScore.Seurat()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScore.md),
one per supported value of its `method` argument. They share the same
interface and only differ in the backend they wrap, hence the common
documentation.

`BenchGeneScoreSeurat.<method>()` is not meant to be called directly: it
is reached through the class that
[`BenchGeneScore.Seurat()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScore.md)
sets on `method` before dispatching.

## Usage

``` r
BenchGeneScoreSeurat.ssGSEA(
  obj,
  pathway,
  assay = NULL,
  layer = "data",
  ...,
  verbose = getFuncOption("verbose")
)

BenchGeneScoreSeurat.AUCell(
  obj,
  pathway,
  assay = NULL,
  layer = "data",
  ...,
  verbose = getFuncOption("verbose")
)

BenchGeneScoreSeurat.UCell(
  obj,
  pathway,
  assay = NULL,
  layer = "data",
  ...,
  verbose = getFuncOption("verbose")
)

BenchGeneScoreSeurat.AddModuleScore(
  obj,
  pathway,
  assay = NULL,
  layer = "data",
  ...,
  verbose = getFuncOption("verbose")
)
```

## Arguments

- obj:

  A `Seurat` object.

- ...:

  Additional arguments forwarded to the wrapped backend.

- verbose:

  Logical. Whether to print progress information. Default:
  `getFuncOption("verbose")`.

## Value

The score matrix of the wrapped backend, with gene sets in rows and
cells in columns, except `"AUCell"`, which returns the
[`AUCell::AUCell_calcAUC()`](https://rdrr.io/pkg/AUCell/man/AUCell_calcAUC.html)
object whose AUC matrix uses the same layout.

## Details

Expression values are always taken with
[`SeuratObject::LayerData()`](https://satijalab.github.io/seurat-object/reference/Layers.html),
so the backend receives a plain gene-by-cell matrix. Each backend then
does its own housekeeping:

- `ssGSEA` splits `...` between
  [`GSVA::ssgseaParam()`](https://rdrr.io/pkg/GSVA/man/ssgseaParam-class.html)
  and [`GSVA::gsva()`](https://rdrr.io/pkg/GSVA/man/gsva.html);

- `AUCell` splits `...` between
  [`AUCell::AUCell_buildRankings()`](https://rdrr.io/pkg/AUCell/man/AUCell_buildRankings.html)
  and
  [`AUCell::AUCell_calcAUC()`](https://rdrr.io/pkg/AUCell/man/AUCell_calcAUC.html),
  and builds the rankings with `plotStats = FALSE`;

- `UCell` passes everything to
  [`UCell::ScoreSignatures_UCell()`](https://rdrr.io/pkg/UCell/man/ScoreSignatures_UCell.html)
  with `name = ""` so that gene-set names stay untouched, then
  transposes the cell-by-gene-set result;

- `AddModuleScore` keeps only the genes present in `assay`, drops the
  gene sets that become empty, and aborts when none of them is left.

## See also

[`BenchGeneScore()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScore.md)
