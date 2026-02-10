# Filter tumor cells (internal)

An internal function that filters tumor cells from a Seurat object based
on metadata column values. This function identifies tumor cells using
pattern matching on cell type labels and creates a subset containing
only tumor cells. It also records dimension information before and after
filtering for traceability.

## Usage

``` r
FilterTumorCell(obj, column2only_tumor = NULL, verbose = TRUE)
```

## Arguments

- obj:

  Seurat object with a column to filter out tumor cells.

- column2only_tumor:

  Name of the column to filter out tumor cells.

- verbose:

  logical, whether to print progress messages

## Value

A Seurat object containing only tumor cells, with the following
attributes stored in `@misc`:

- `self_dim`: Dimensions of the filtered object

- `raw_dim`: Original dimensions before filtering

- `column2only_tumor`: The column name used for filtering

If `column2only_tumor` is `NULL` or the specified column is not found,
returns the original object unchanged.

## See also

Other single_cell_preprocess:
[`FindRobustElbow()`](https://wanglabcsu.github.io/sigbridger/reference/FindRobustElbow.md),
[`Pattern2Colname()`](https://wanglabcsu.github.io/sigbridger/reference/Pattern2Colname.md),
[`QCPatternDetect()`](https://wanglabcsu.github.io/sigbridger/reference/QCPatternDetect.md),
[`RegisterSeuratMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterSeuratMethod.md),
[`SCAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/SCAnnotate.md),
[`SCIntegrate()`](https://wanglabcsu.github.io/sigbridger/reference/SCIntegrate.md),
[`SCPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcess.md),
[`SCPreProcessStrategy`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcessStrategy.md),
[`compatible_with_3.0.2()`](https://wanglabcsu.github.io/sigbridger/reference/compatible_with_3.0.2.md)
