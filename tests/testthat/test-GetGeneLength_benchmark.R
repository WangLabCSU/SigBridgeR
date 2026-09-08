skip("skipped GetGeneLength benchmarking test")
skip_if_not_installed("GenomicRanges")
skip_if_not_installed("rtracklayer")


test_that("Compare computation accuracy", {
  gtf_path <- "/data/resource/wanglab/GRCh37/Homo_sapiens.GRCh37.75.gtf.gz"

  gtf <- rtracklayer::import(gtf_path)
  exon <- gtf[gtf$type == "exon"]
  exon$gene_id <- sub("\\..*$", "", exon$gene_id)
  gene_length <- tapply(seq_along(exon), exon$gene_id, function(i) {
    sum(GenomicRanges::width(GenomicRanges::reduce(exon[i])))
  })

  gene_length <- data.frame(
    gene_id = names(gene_length),
    gene_length = as.numeric(gene_length),
    row.names = NULL
  )

  qs2::qd_save(
    gene_length,
    "/data/resource/wanglab/GRCh37/gene_length_r_computed.qs2"
  )

  gene_length_cpp <- GetGeneLength(gtf_path)
  gene_length_cpp <- gene_length_cpp[order(names(gene_length_cpp))]
  gene_length_cpp <- data.frame(
    gene_id = names(gene_length_cpp),
    gene_length = as.numeric(gene_length_cpp),
    row.names = NULL
  )

  qs2::qd_save(
    gene_length_cpp,
    "/data/resource/wanglab/GRCh37/gene_length_cpp_computed.qs2"
  )

  expect_lt(max(gene_length$gene_length - gene_length_cpp$gene_length), 1e-6)
})
