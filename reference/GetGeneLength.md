# Compute Gene Lengths from a GTF File or GTF Lines

An S3 generic that computes per-gene lengths by parsing `exon` features
from either a GTF/GTF.GZ file or a character vector of GTF lines. The
reported length for each gene is the total merged (non-overlapping) span
of its exons.

The `character` method accepts two forms of input:

- A **single-element** character vector `x` is treated as a GTF file
  path (optionally gzip-compressed). The file is streamed line-by-line
  and is never loaded fully into memory.

- A **multi-element** character vector `x` is treated as the GTF lines
  themselves (e.g. the output of
  [`readLines()`](https://rdrr.io/r/base/readLines.html)).

Overlapping or adjacent exon intervals are merged within each
(`gene_id`, chromosome, strand) group, and the length of the merged span
is reported. Ensembl `gene_id` version suffixes (e.g.
`ENSG00000139618.15`) are preserved, so gene identifiers from different
annotation builds remain distinct.

By default, exon records on **every** chromosome are kept. Set
`canonical_chr_only = TRUE` to restrict the analysis to canonical
chromosomes (autosomes `1-22` and chromosomes `X`, `Y`, and `MT`) and
ignore any other contigs or scaffolds present in the annotation.

## Usage

``` r
GetGeneLength(x, verbose = TRUE, canonical_chr_only = FALSE, ...)

# Default S3 method
GetGeneLength(x, verbose = TRUE, canonical_chr_only = FALSE, ...)

# S3 method for class 'character'
GetGeneLength(x, verbose = TRUE, canonical_chr_only = FALSE, ...)
```

## Arguments

- x:

  The object to dispatch on:

  - A length-1 character string: path to a GTF or GTF.GZ file.

  - A character vector with more than one element: the GTF lines.

- verbose:

  A logical flag. If `TRUE` (default), cli-style progress messages (with
  ANSI colors) are printed while reading, parsing, and merging the exon
  intervals.

- canonical_chr_only:

  A logical flag. If `TRUE`, only `exon` records located on canonical
  chromosomes (autosomes `1-22`, and chromosomes `X`, `Y`, `MT`) are
  retained; exons on any other contig/scaffold are ignored. If `FALSE`
  (default), exons on every chromosome are considered.

- ...:

  Reserved for future use.

## Value

A named numeric vector of gene lengths, one entry per gene, with the
gene ids as names. An empty (zero-length) named vector is returned when
the input contains no valid `exon` records.

## Examples

``` r
if (FALSE) { # \dontrun{
# Gene lengths from a gzip-compressed GTF file
lengths <- GetGeneLength("Homo_sapiens.GRCh37.75.gtf.gz")

# Gene lengths from GTF lines already read into R
x <- readLines("Homo_sapiens.GRCh37.75.gtf.gz")
lengths <- GetGeneLength(x)
} # }
```
