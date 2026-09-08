#' Compute Gene Lengths from a GTF File or GTF Lines
#'
#' @description
#' An S3 generic that computes per-gene lengths by parsing `exon` features
#' from either a GTF/GTF.GZ file or a character vector of GTF lines. The
#' reported length for each gene is the total merged (non-overlapping) span
#' of its exons.
#'
#' The `character` method accepts two forms of input:
#' * A **single-element** character vector `x` is treated as a GTF file path
#'   (optionally gzip-compressed). The file is streamed line-by-line and is
#'   never loaded fully into memory.
#' * A **multi-element** character vector `x` is treated as the GTF lines
#'   themselves (e.g. the output of [readLines()]).
#'
#' Overlapping or adjacent exon intervals are merged within each (`gene_id`,
#' chromosome, strand) group, and the length of the merged span is reported.
#' Ensembl `gene_id` version suffixes (e.g. `ENSG00000139618.15`) are preserved,
#' so gene identifiers from different annotation builds remain distinct.
#'
#' By default, exon records on **every** chromosome are kept. Set
#' `canonical_chr_only = TRUE` to restrict the analysis to canonical
#' chromosomes (autosomes `1-22` and chromosomes `X`, `Y`, and `MT`) and ignore
#' any other contigs or scaffolds present in the annotation.
#'
#' @param x The object to dispatch on:
#'   * A length-1 character string: path to a GTF or GTF.GZ file.
#'   * A character vector with more than one element: the GTF lines.
#' @param verbose A logical flag. If `TRUE` (default), cli-style progress
#'   messages (with ANSI colors) are printed while reading, parsing, and
#'   merging the exon intervals.
#' @param canonical_chr_only A logical flag. If `TRUE`, only `exon` records
#'   located on canonical chromosomes (autosomes `1-22`, and chromosomes `X`,
#'   `Y`, `MT`) are retained; exons on any other contig/scaffold are ignored.
#'   If `FALSE` (default), exons on every chromosome are considered.
#' @param ... Reserved for future use.
#'
#' @return A named numeric vector of gene lengths, one entry per gene, with
#'   the gene ids as names. An empty (zero-length) named vector is returned
#'   when the input contains no valid `exon` records.
#'
#' @examples
#' \dontrun{
#' # Gene lengths from a gzip-compressed GTF file
#' lengths <- GetGeneLength("Homo_sapiens.GRCh37.75.gtf.gz")
#'
#' # Gene lengths from GTF lines already read into R
#' x <- readLines("Homo_sapiens.GRCh37.75.gtf.gz")
#' lengths <- GetGeneLength(x)
#' }
#'
#' @export
GetGeneLength <- function(x, verbose = TRUE, canonical_chr_only = FALSE, ...) {
  chk::chk_flag(verbose)
  chk::chk_flag(canonical_chr_only)
  UseMethod("GetGeneLength")
}


#' @rdname GetGeneLength
#' @export
GetGeneLength.default <- function(
  x,
  verbose = TRUE,
  canonical_chr_only = FALSE,
  ...
) {
  available_methods <- gsub(
    "^GetGeneLength\\.",
    "",
    utils::methods(GetGeneLength)
  )
  Abort(
    "GetGeneLength() is not implemented for objects of {.cls {class(x)}}",
    "Expected class: {.cls {available_methods}}",
    type = "[CLASS ERROR]"
  )
}

#' @rdname GetGeneLength
#' @export
GetGeneLength.character <- function(
  x,
  verbose = TRUE,
  canonical_chr_only = FALSE,
  ...
) {
  if (length(x) == 1L) {
    return(GetGeneLength_file_impl(
      x = x,
      verbose = verbose,
      canonical_chr_only = canonical_chr_only,
      ...
    ))
  }
  GetGeneLength_lines_impl(
    x = x,
    verbose = verbose,
    canonical_chr_only = canonical_chr_only,
    ...
  )
}

GetGeneLength_file_impl <- function(
  x,
  verbose = TRUE,
  canonical_chr_only = FALSE,
  ...
) {
  check_dots_empty()
  if (!file.exists(x)) {
    Abort("{.file {x}} is not found", type = "[FILE ERROR]")
  } else if (!endsWith(x, ".gz")) {
    Abort(
      "{.file {x}} is not a gzipped file",
      "E.g. {.path Homo_sapiens.GRCh37.75.gtf.gz}",
      type = "[FILE FORMAT]"
    )
  }
  gtf_file_to_gene_length(
    path = x,
    verbose = verbose,
    canonical_chr_only = canonical_chr_only
  )
}


GetGeneLength_lines_impl <- function(
  x,
  verbose = TRUE,
  canonical_chr_only = FALSE,
  ...
) {
  check_dots_empty()
  r_gtf_lines_to_gene_length(
    lines = x,
    verbose = verbose,
    canonical_chr_only = canonical_chr_only
  )
}
