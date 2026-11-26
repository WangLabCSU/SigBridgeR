# Single-sample Gene Set Scoring for Single-Cell Data

Scores one or more gene sets in every cell of a screened single-cell
dataset with `ssGSEA`, `AUCell`, `UCell` or `AddModuleScore`. All
backends return the same gene-set-by-cell matrix layout, so the
cell-level signal recovered by different screening methods can be
compared on a common ground.

## Usage

``` r
BenchGeneScore(obj, ...)

# S3 method for class '`SigBridgeR::BenchGeneScoreConfig`'
BenchGeneScore(obj, ..., verbose = getFuncOption("verbose"))

# S3 method for class '`SigBridgeR::ScreenMethodResult`'
BenchGeneScore(obj, ...)

# S3 method for class 'Seurat'
BenchGeneScore(
  obj,
  method = c("ssGSEA", "AUCell", "UCell", "AddModuleScore"),
  ...
)
```

## Arguments

- obj:

  The object to score, one of:

  - a `Seurat` object, scored as is;

  - a
    [ScreenMethodResult](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodResult.md)
    object, whose `scRNA_data` slot is scored;

  - a
    [BenchGeneScoreConfig](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfig.md)
    object, whose `screen_method` entries are all run through
    [`Screen()`](https://wanglabcsu.github.io/sigbridger/reference/Screen.md)
    first and the accumulated labels are scored in a single pass.

- ...:

  Additional named arguments forwarded to the selected backend.
  Arguments are matched against the formal argument names of the backend
  and split automatically over the two stages when a method needs two
  calls:

  - `"ssGSEA"`:
    [`GSVA::ssgseaParam()`](https://rdrr.io/pkg/GSVA/man/ssgseaParam-class.html)
    and [`GSVA::gsva()`](https://rdrr.io/pkg/GSVA/man/gsva.html);

  - `"AUCell"`:
    [`AUCell::AUCell_buildRankings()`](https://rdrr.io/pkg/AUCell/man/AUCell_buildRankings.html)
    and
    [`AUCell::AUCell_calcAUC()`](https://rdrr.io/pkg/AUCell/man/AUCell_calcAUC.html);

  - `"UCell"`:
    [`UCell::ScoreSignatures_UCell()`](https://rdrr.io/pkg/UCell/man/ScoreSignatures_UCell.html);

  - `"AddModuleScore"`:
    [`Seurat::AddModuleScore()`](https://satijalab.org/seurat/reference/AddModuleScore.html).

- verbose:

  Logical. Whether to report progress information. Default: `TRUE`.

- method:

  Character. Scoring backend, one of `"ssGSEA"`, `"AUCell"`, `"UCell"`
  or `"AddModuleScore"`. Default: `"ssGSEA"`.

## Value

A numeric matrix with gene sets in rows and cells in columns. Gene sets
dropped by backend-specific filtering are omitted from the result.

## Details

### Dispatch

`BenchGeneScore()` dispatches on the class of `obj`, so a
[BenchGeneScoreConfig](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfig.md)
is benchmarked and a
[ScreenMethodResult](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodResult.md)
is unwrapped before the actual scoring happens.
`BenchGeneScore.Seurat()` is itself a generic dispatching on `method`:
every backend lives in its own `BenchGeneScoreSeurat.<method>` function
and its package is loaded on demand by
[`rlang::check_installed()`](https://rlang.r-lib.org/reference/is_installed.html),
so neither `GSVA`, `AUCell`, `UCell` nor `Seurat` is required to merely
call the generic.

### Config input

For a
[BenchGeneScoreConfig](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfig.md)
input, the configured `screen_method` entries are run one after another
with
[`Screen()`](https://wanglabcsu.github.io/sigbridger/reference/Screen.md).
Every run adds its labels to the `meta.data` of the working `Seurat`
object and the labels are accumulated across methods, after which the
scores are computed once. Scores are therefore comparable across
screening methods by construction.

### Backend notes

- `"UCell"` is called with `name = ""` so that gene-set names stay
  untouched, and its cell-by-gene-set output is transposed to the common
  layout.

- `"AddModuleScore"` drops the genes absent from `assay`, drops the gene
  sets that become empty, and stores the scores under an internal
  `.BenchGeneScore` prefix so that the `meta.data` of `obj` is left
  untouched. It aborts when no gene set survives the filtering.

- `"AUCell"` builds the rankings with `plotStats = FALSE`, since the
  plot is only meaningful for a manual inspection.

## See also

[`BenchGeneScoreConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfig.md),
[`BenchGeneScoreConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfigSet.md),
[`Screen()`](https://wanglabcsu.github.io/sigbridger/reference/Screen.md),
[`BenchGSEA()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEA.md),
[`BenchAcc()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAcc.md)

Other Benchmark:
[`BenchAcc()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAcc.md),
[`BenchGSEA()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEA.md),
[`BenchMemory()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemory.md)

## Examples

``` r
# \donttest{
counts <- matrix(
  rep(seq_len(20), times = 5) + rep(seq_len(5), each = 20),
  nrow = 20,
  dimnames = list(paste0("Gene", seq_len(20)), paste0("Cell", seq_len(5)))
)
obj <- Seurat::CreateSeuratObject(counts = counts)
#> Warning: Data is of class matrix. Coercing to dgCMatrix.
obj <- Seurat::NormalizeData(obj, verbose = FALSE)
pathways <- list(
  First = paste0("Gene", 1:5),
  Second = paste0("Gene", 6:10)
)
scores <- try(BenchGeneScore(obj, pathways, method = "ssGSEA", verbose = FALSE)) # require GSVA
#> ℹ 100 nonzeros (less than 2^31) and 0.00% sparsity
#> ! Input expression data is sparse, but the ssGSEA algorithm does not deal with sparsity in any specific way, and data will be converted into a dense matrix format
# }
```
