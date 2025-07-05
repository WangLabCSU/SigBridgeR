library(dplyr)
setwd("~/R/Project/R_code/SigBridgeR")
source("./R/11-PreProcess.R")
source("./R/90-Test.R")

tcga_exp_count <- readRDS(
  "~/R/Project/R_code/sideWORK/Scissored-GC-ecDNA/data/input/TCGA-LUAD.exp.count.rds"
) %>%
  magrittr::set_colnames(substr(colnames(.), 1, 15)) %>%
  BulkPreProcess()

tcga_ms_indel = read.csv(
  "~/R/Project/R_code/WORK/data/input_data/TCGA_WES_sigProfiler_ID_signatures_in_samples.csv"
) %>%
  MSPreProcess(
    col_thresh = 0.05
  )

GSE150290 = readRDS(
  "~/Project/Scissored-GC-ecDNA/data/input/GSE150290.rds"
) |>
  Seurat::FindNeighbors(dims = 1:16) |>
  Seurat::FindClusters(resolution = 0.7) |>
  Seurat::RunTSNE(dims = 1:16) |>
  Seurat::RunUMAP(dims = 1:16)

col_id_seq = which(IdentifyDataColumn(tcga_ms_sbs))

col_id = 3

match_result = MatchSample(
  ms_signature = tcga_ms_sbs,
  TCGA_exp_count = tcga_exp_count,
  col_id = col_id
)

scissor_result = DoScissor(
  matched_bulk = match_result$matched_TCGA_exp_count,
  sc_data = GSE150290,
  phenotype = match_result$phenotype,
  label_type = match_result$ms_select,
  dir2save_scissor_inputs = "test_result",
  scissor_family = "binomial",
  reliability_test = FALSE,
  nfold = 5
)
scissor_result = DoScissor(
  path2load_scissor_cache = "test_result/Scissor_inputs.RData",
  sc_data = GSE150290,
  label_type = match_result$ms_select,
  scissor_family = "binomial",
  reliability_test = FALSE,
  nfold = 5
)
