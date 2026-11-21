# Tests for src/GetGeneLength.cpp — gtf_file_to_gene_length() and
# r_gtf_lines_to_gene_length(), exercised through their R (Rcpp-exported)
# interfaces and the public GetGeneLength() S3 generic.
#
# Run independently with:
#   devtools::load_all()
#   testthat::test_file("tests/testthat/test-GetGeneLength_cpp.R")
#
# All GTF inputs are tiny (1-6 lines) to keep the suite fast.

# Helper: build a single GTF line with 9 tab-separated fields.
.gtf_line <- function(chr = "1", feature = "exon", start, end, strand = "+", gene_id = "G", extra = "") {
  paste(
    chr,
    "test_src",
    feature,
    start,
    end,
    ".",
    strand,
    ".",
    paste0('gene_id "', gene_id, '";', extra),
    sep = "\t"
  )
}

# ---- r_gtf_lines_to_gene_length(): merging ----------------------------------

describe("r_gtf_lines_to_gene_length - exon merging", {
  it("reports the length of a single exon", {
    lines <- .gtf_line(start = 100, end = 200, gene_id = "G1")
    res <- suppressMessages(r_gtf_lines_to_gene_length(lines, verbose = FALSE))
    expect_equal(res, c(G1 = 101))
  })

  it("merges overlapping exons of the same gene", {
    lines <- c(
      .gtf_line(start = 100, end = 200, gene_id = "G1"),
      .gtf_line(start = 150, end = 250, gene_id = "G1")
    )
    res <- suppressMessages(r_gtf_lines_to_gene_length(lines, verbose = FALSE))
    expect_equal(res, c(G1 = 151))
  })

  it("sums non-overlapping exons of the same gene", {
    lines <- c(
      .gtf_line(start = 100, end = 200, gene_id = "G1"),
      .gtf_line(start = 300, end = 400, gene_id = "G1")
    )
    res <- suppressMessages(r_gtf_lines_to_gene_length(lines, verbose = FALSE))
    expect_equal(res, c(G1 = 202))
  })

  it("computes multiple genes independently and sorts by gene_id", {
    lines <- c(
      .gtf_line(start = 20, end = 29, gene_id = "Z"),
      .gtf_line(start = 1, end = 10, gene_id = "A"),
      .gtf_line(start = 50, end = 55, gene_id = "M")
    )
    res <- suppressMessages(r_gtf_lines_to_gene_length(lines, verbose = FALSE))
    expect_equal(res, c(A = 10, M = 6, Z = 10))
    expect_identical(names(res), c("A", "M", "Z"))
  })

  it("preserves gene_id version suffixes", {
    lines <- .gtf_line(start = 1, end = 10, gene_id = "ENSG00000139618.15")
    res <- suppressMessages(r_gtf_lines_to_gene_length(lines, verbose = FALSE))
    expect_equal(res, c(`ENSG00000139618.15` = 10))
  })

  it("accepts unquoted gene_id values", {
    lines <- sub(
      'gene_id "G1";',
      "gene_id G1;",
      .gtf_line(start = 1, end = 10, gene_id = "G1"),
      fixed = TRUE
    )
    res <- suppressMessages(r_gtf_lines_to_gene_length(lines, verbose = FALSE))
    expect_equal(res, c(G1 = 10))
  })

  it("sums lengths across chromosomes and strands for the same gene", {
    lines <- c(
      .gtf_line(chr = "1", strand = "+", start = 1, end = 10, gene_id = "G1"),
      .gtf_line(chr = "2", strand = "-", start = 1, end = 10, gene_id = "G1")
    )
    res <- suppressMessages(r_gtf_lines_to_gene_length(lines, verbose = FALSE))
    expect_equal(res, c(G1 = 20))
  })
})

# ---- r_gtf_lines_to_gene_length(): filtering / edge cases -------------------

describe("r_gtf_lines_to_gene_length - filtering and edge cases", {
  it("skips comment lines and non-exon features", {
    lines <- c(
      "# a comment",
      .gtf_line(feature = "gene", start = 1, end = 100, gene_id = "G1"),
      .gtf_line(feature = "exon", start = 10, end = 20, gene_id = "G1")
    )
    res <- suppressMessages(r_gtf_lines_to_gene_length(lines, verbose = FALSE))
    expect_equal(res, c(G1 = 11))
  })

  it("filters non-canonical contigs when canonical_chr_only = TRUE", {
    lines <- c(
      .gtf_line(chr = "1", start = 1, end = 10, gene_id = "A"),
      .gtf_line(chr = "scaffold_123", start = 1, end = 10, gene_id = "B")
    )
    res <- suppressMessages(
      suppressWarnings(
        r_gtf_lines_to_gene_length(lines, verbose = FALSE, canonical_chr_only = TRUE)
      )
    )
    expect_equal(res, c(A = 10))
  })

  it("accepts chr-prefixed and X/Y/M/MT chromosomes in canonical mode", {
    lines <- c(
      .gtf_line(chr = "chr1", start = 1, end = 10, gene_id = "A"),
      .gtf_line(chr = "chrX", start = 1, end = 10, gene_id = "B"),
      .gtf_line(chr = "MT", start = 1, end = 10, gene_id = "C")
    )
    res <- suppressMessages(
      r_gtf_lines_to_gene_length(lines, verbose = FALSE, canonical_chr_only = TRUE)
    )
    expect_setequal(names(res), c("A", "B", "C"))
  })

  it("rejects leading-zero chromosome names in canonical mode", {
    lines <- .gtf_line(chr = "01", start = 1, end = 10, gene_id = "G1")
    expect_warning(
      res <- suppressMessages(
        r_gtf_lines_to_gene_length(lines, verbose = FALSE, canonical_chr_only = TRUE)
      ),
      "No exon with a valid gene_id"
    )
    expect_length(res, 0)
  })

  it("skips NA lines", {
    lines <- c(NA_character_, .gtf_line(start = 1, end = 10, gene_id = "G1"))
    res <- suppressMessages(r_gtf_lines_to_gene_length(lines, verbose = FALSE))
    expect_equal(res, c(G1 = 10))
  })

  it("warns about malformed coordinates and skips those exons", {
    lines <- c(
      .gtf_line(start = 1, end = 10, gene_id = "G1"),
      .gtf_line(start = "abc", end = 20, gene_id = "G2")
    )
    expect_warning(
      res <- suppressMessages(r_gtf_lines_to_gene_length(lines, verbose = FALSE)),
      "1 exon lines with invalid coordinates"
    )
    expect_equal(res, c(G1 = 10))
  })

  it("warns and returns an empty named vector when no valid exon exists", {
    lines <- .gtf_line(feature = "gene", start = 1, end = 10, gene_id = "G1")
    expect_warning(
      res <- suppressMessages(r_gtf_lines_to_gene_length(lines, verbose = FALSE)),
      "No exon with a valid gene_id"
    )
    expect_length(res, 0)
    expect_type(res, "double")
  })

  it("emits cli progress messages when verbose = TRUE", {
    lines <- .gtf_line(start = 1, end = 10, gene_id = "G1")
    expect_message(
      res <- r_gtf_lines_to_gene_length(lines, verbose = TRUE),
      "Computed merged gene lengths"
    )
    expect_equal(res, c(G1 = 10))
  })
})

# ---- gtf_file_to_gene_length(): file input ----------------------------------

describe("gtf_file_to_gene_length - file input", {
  it("reads a plain-text GTF file", {
    path <- withr::local_tempfile(fileext = ".gtf")
    writeLines(.gtf_line(start = 1, end = 10, gene_id = "G1"), path)
    res <- suppressMessages(gtf_file_to_gene_length(path, verbose = FALSE))
    expect_equal(res, c(G1 = 10))
  })

  it("reads a gzipped GTF file", {
    path <- withr::local_tempfile(fileext = ".gtf.gz")
    con <- gzfile(path, "w")
    writeLines(.gtf_line(start = 1, end = 10, gene_id = "G1"), con)
    close(con)
    res <- suppressMessages(gtf_file_to_gene_length(path, verbose = FALSE))
    expect_equal(res, c(G1 = 10))
  })

  it("merges exons across a multi-line file", {
    path <- withr::local_tempfile(fileext = ".gtf")
    writeLines(
      c(
        .gtf_line(start = 100, end = 200, gene_id = "G1"),
        .gtf_line(start = 150, end = 250, gene_id = "G1"),
        .gtf_line(start = 1, end = 5, gene_id = "G2")
      ),
      path
    )
    res <- suppressMessages(gtf_file_to_gene_length(path, verbose = FALSE))
    expect_equal(res, c(G1 = 151, G2 = 5))
  })

  it("errors on a missing file", {
    expect_error(
      gtf_file_to_gene_length("/nonexistent/path.gtf", verbose = FALSE),
      "Cannot open file"
    )
  })
})

# ---- GetGeneLength(): public S3 interface -----------------------------------

describe("GetGeneLength - public S3 interface", {
  it("accepts a character vector of GTF lines", {
    lines <- c(
      .gtf_line(start = 1, end = 10, gene_id = "A"),
      .gtf_line(start = 20, end = 25, gene_id = "B")
    )
    res <- suppressMessages(GetGeneLength(lines, verbose = FALSE))
    expect_equal(res, c(A = 10, B = 6))
  })

  it("accepts a gzipped GTF file path", {
    path <- withr::local_tempfile(fileext = ".gtf.gz")
    con <- gzfile(path, "w")
    writeLines(.gtf_line(start = 1, end = 10, gene_id = "G1"), con)
    close(con)
    res <- suppressMessages(GetGeneLength(path, verbose = FALSE))
    expect_equal(res, c(G1 = 10))
  })

  it("rejects a non-gzipped file path", {
    path <- withr::local_tempfile(fileext = ".gtf")
    writeLines(.gtf_line(start = 1, end = 10, gene_id = "G1"), path)
    expect_error(GetGeneLength(path, verbose = FALSE), "not a gzipped file")
  })

  it("rejects objects with no method", {
    expect_error(GetGeneLength(1:3), "not implemented")
  })
})
