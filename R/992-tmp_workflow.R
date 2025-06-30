library(dplyr)

tcga_exp_count <- readRDS(
  "~/Data/TCGA-LUAD.exp.count.rds"
) %>%
  magrittr::set_colnames(substr(colnames(.), 1, 15)) %>%
  BulkPreProcess()

tcga_ms_sbs = read.csv(
  "~/Data/MutationalSignature/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv"
) %>%
  MSPreProcess(
    thresh = 0.05
  )


