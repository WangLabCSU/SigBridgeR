# Safely Add Miscellaneous Data to Seurat Object

Adds arbitrary data to the `@misc` slot of a Seurat object with
automatic key conflict resolution. If the key already exists,
automatically appends a numeric suffix to ensure unique key naming
(e.g., "mykey_1", "mykey_2").

## Usage

``` r
AddMisc(seurat_obj, ..., cover = TRUE)
```

## Arguments

- seurat_obj:

  A Seurat object to modify

- ...:

  key-value pairs to add to the `@misc` slot.

- cover:

  Logical indicating whether to overwrite existing data. If (default
  TRUE).

## Value

The modified Seurat object with added `@misc` data. The original object
structure is preserved with no other modifications.

## Key Generation Rules

1.  If `key` doesn't exist: uses as-is

2.  If `key` exists: appends the next available number (e.g., "key_1",
    "key_2")

3.  If numbered keys exist (e.g., "key_2"): increments the highest
    number

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
seurat_obj <- AddMisc(seurat_obj, "QC_stats" = qc_df)

# Auto-incrementing example
seurat_obj <- AddMisc(seurat_obj, markers = markers1)
seurat_obj <- AddMisc(seurat_obj, markers = markers2, cover=FALSE)
# Stores as "markers" and "markers_1"

} # }
```
