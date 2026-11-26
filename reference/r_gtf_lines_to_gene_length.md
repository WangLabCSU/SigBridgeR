# Compute Gene Lengths from GTF Lines

Applies the same exon-merging logic as
[`gtf_file_to_gene_length()`](https://wanglabcsu.github.io/sigbridger/reference/gtf_file_to_gene_length.md)
to GTF lines already held in memory.

## Usage

``` r
r_gtf_lines_to_gene_length(lines, verbose = TRUE, canonical_chr_only = FALSE)
```

## Arguments

- lines:

  Character vector of GTF lines.

- verbose:

  `logical`. Whether to emit progress messages. Default: `TRUE`.

- canonical_chr_only:

  `logical`. Whether to keep only canonical chromosomes. Default:
  `FALSE`.

## Value

A numeric vector named by `gene_id` and sorted by it, holding merged
exon lengths.
