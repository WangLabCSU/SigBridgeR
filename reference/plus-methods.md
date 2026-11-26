# Merge two objects with `+`

`x + y` combines two screening results into a single Seurat object via
[`MergeResult()`](https://wanglabcsu.github.io/sigbridger/reference/MergeResult.md).
Both operands may be
[ScreenMethodResult](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodResult.md)
objects or Seurat objects; a
[ScreenMethodResult](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodResult.md)
is unwrapped to its `scRNA_data` slot before merging.

## Usage

``` r
# S3 method for class '`SigBridgeR::ScreenMethodResult`'
x + y

# S3 method for class 'Seurat'
x + y
```

## Arguments

- x:

  A ScreenMethodResult or Seurat object as the main receiver.

- y:

  A ScreenMethodResult or Seurat object as the main donor.

## Value

A merged Seurat object.
