# Compute Gene Lengths from a GTF File

Streams a plain or gzipped GTF file and merges overlapping exons by
gene.

## Usage

``` r
gtf_file_to_gene_length(path, verbose = TRUE, canonical_chr_only = FALSE)
```

## Arguments

- path:

  `character`. Path to the GTF file, optionally `.gz`.

- verbose:

  `logical`. Whether to emit progress messages. Default: `TRUE`.

- canonical_chr_only:

  `logical`. Whether to keep only canonical chromosomes. Default:
  `FALSE`.

## Value

A numeric vector named by `gene_id` and sorted by it, holding merged
exon lengths.
