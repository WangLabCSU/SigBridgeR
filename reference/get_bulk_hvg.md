# Select Highly Variable Genes from Bulk Data

Ranks genes by the standard deviation of their row in `patdata`.

## Usage

``` r
get_bulk_hvg(patdata, common_genes, n_hvg)
```

## Arguments

- patdata:

  Numeric matrix with genes in rows.

- common_genes:

  Character vector of gene names, one per row of `patdata`.

- n_hvg:

  `integer`. Maximum number of genes to return.

## Value

Character vector of at most `n_hvg` gene names, ordered by decreasing
row standard deviation; `character(0)` when there is no valid input.
