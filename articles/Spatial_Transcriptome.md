# Application of SigBridgeR to Single-cell Spatial Transcriptomics Data

``` r

library(SigBridgeR)
#> ✔ SigBridgeR v3.7.0 loaded
library(Seurat)
#> Loading required package: SeuratObject
#> Loading required package: sp
#> 
#> Attaching package: 'SeuratObject'
#> The following objects are masked from 'package:base':
#> 
#>     intersect, t
```

## Load Data

### Single Cell Transcriptomics

We will use `GSE274103` from NCBI GEO as an example. The file
organization process is not described in detail here.

``` r

data_dir <- "GSE274103"

samples <- c(
  "GSM8443449_PDAC-p1",
  "GSM8443450_PDAC-p2",
  "GSM8443451_PDAC-p3",
  "GSM8443452_PDAC-p4",
  "GSM8443453_PDAC-p5"
)

spatial_list <- lapply(samples, function(sample) {
  seurat <- Load10X_Spatial(
    data.dir = file.path(data_dir, sample),
    slice = sample,
    assay = "Spatial"
  )

  seurat <- subset(
    seurat,
    nCount_Spatial > 0L &
      nFeature_Spatial > 0L
  )

  SCTransform(seurat, assay = "Spatial")
})

spatials <- merge(
  x = spatial_list[[1]],
  y = unlist(spatial_list[-1]),
  add.cell.ids = samples,
  merge.data = TRUE,
  project = "Integrated Spatial Seurat"
)

# DefaultAssay(spatials) <- "SCT"
VariableFeatures(spatials) <- lapply(spatial_list, function(seurat) {
  VariableFeatures(seurat)
}) %>%
  unlist() %>%
  unique()
spatials <- RunPCA(spatials) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30)

# SeuratObject::SaveSeuratRds(spatials, "seurat.rds")
```

In theory, this workflow could be simplified using `SCPreProcess`;
however, due to the performance issue reported in [Seurat
\#10153](https://github.com/satijalab/seurat/issues/10153)—where
`SCTransform` hangs or slows down significantly when called via
`do.call` (as `SCPreProcess` does internally)—we instead use a custom
workflow here to avoid the slowdown.

(It hasn’t been fixed yet. If it gets resolved, please kindly notify me
in the issue—thank you! :))

### Bulk RNA Seq & Survival

The bulk RNA-seq example data used here is from TCGA-PAAD. Due to
copyright restrictions, it is provided solely for illustrative purposes
in the vignette.

``` r

bulk <- readRDS("TCGA-PAAD.rds")
#        TCGA-2J-AAB1-01 TCGA-2J-AAB4-01 TCGA-2J-AAB6-01 TCGA-2J-AAB8-01
# TSPAN6       10.336507       10.992938       10.143383        9.415742
# TNMD          0.000000        0.000000        0.000000        2.584963
# DPM1          9.967226       10.327553       10.974415        9.527477
# SCYL3         9.479780        9.481799        8.370687        9.142107
```

The corresponding TCGA survival data can be obtained from **UCSC Xena
(via UCSCXenaShiny)**. And we match the corresponding samples in both
the bulk RNA-seq data and the survival information.

``` r

library(UCSCXenaShiny)
tcga_surv <- load_data("tcga_surv")

paad_samples <- colnames(bulk)

tcga_paad_surv <- dplyr::filter(tcga_surv, sample %in% paad_samples) %>%
  dplyr::select(sample, OS.time, OS) %>%
  dplyr::rename(time = OS.time, status = OS) %>%
  dplyr::filter(status != "NA", time != "NA", !is.na(status), !is.na(time)) %>%
  tibble::column_to_rownames(var = "sample")

head(tcga_paad_surv)
#                 time status
# TCGA-2J-AAB1-01   66      1
# TCGA-2J-AAB4-01  729      0
# TCGA-2J-AAB6-01  293      1
# TCGA-2J-AAB8-01   80      0
# TCGA-2J-AAB9-01  627      1
# TCGA-2J-AABA-01  607      1

bulk <- bulk[, rownames(tcga_paad_surv)]
```

We will use a single-cell phenotype-based screening algorithm to
identify cell populations associated with poor survival prognosis.
However, to determine the identity of these cells, we first need to
perform cell type annotation on the data.

## Annotation of Single Cell Data

Here we demonstrate the usage of **mLLMCelltype**. We also provide
support for **SingleR** and **CellTypist** (see [Auxiliary
Utils](https://wanglabcsu.github.io/SigBridgeR/articles/Other_Function_Details.html)
for details and usage).

Here we use DeepSeek v3 only as an example. Generally speaking, the more
powerful the model and the greater the number of models used, the higher
the prediction accuracy of mLLMCelltype.

``` r

spatials <- Seurat::PrepSCTFindMarkers(spatials)

spatials <- SCAnnotate(
  spatials,
  models = c("deepseek-chat"),
  api_keys = list(
    deepseek = "sk-1234567890"
  )
)

# SeuratObject::SaveSeuratRds(spatials, "seurat.rds")
```

``` r

# ℹ [2026/02/09 11:26:14] [mLLMCelltype] Start annotating cell types
# ℹ [2026/02/09 11:26:14] Find marker genes for each clusters
# Calculating cluster 0
# Calculating cluster 1
# Calculating cluster 2
# Calculating cluster 3
# Calculating cluster 4
# Calculating cluster 5
# Calculating cluster 6
# Calculating cluster 7
# Calculating cluster 8
# Calculating cluster 9
# Calculating cluster 10
# Calculating cluster 11
# Calculating cluster 12
# Calculating cluster 13
# Calculating cluster 14
# Calculating cluster 15
# Calculating cluster 16
# Calculating cluster 17
# Calculating cluster 18
# Calculating cluster 19
# Calculating cluster 20
# Calculating cluster 21
# Calculating cluster 22
# Calculating cluster 23
# Calculating cluster 24
# ℹ [2026/02/09 11:28:42] Large language models cell type Annotating
# ###
# # LLM Output
# ###
# ✔ [2026/02/09 11:28:46] Annotation Finished
```

``` r

table(spatials$mllmcelltype_cell_type)
```

``` r

#     Acinar cells           Adipocytes              B cells         Chondrocytes    Endothelial cells
#              594                  633                 2314                  920                 2122
#      Enterocytes     Epithelial cells          Fibroblasts        Gastric cells         Goblet cells
#             2018                 4799                 3395                  547                  646
#    Keratinocytes          Macrophages    Mesothelial cells Neuroendocrine cells  Smooth muscle cells
#             1135                  879                  242                 1751                  254
#          T cells
#             1187
```

## Screen Phenotype-associated Cells

We can now run the screening. Let’s try `Scissor`.

``` r

res <- Screen(
  matched_bulk = bulk,
  sc_data = spatials,
  phenotype = tcga_paad_surv,
  phenotype_class = "survival",
  screen_method = "Scissor",
  assay = "SCT"
)

# SeuratObject::SaveSeuratRds(res$scRNA_data, "screened_seurat.rds")
```

``` r

# ℹ `label_type` not specified or not of length 1, using "Scissor"
# ℹ [2026/02/09 15:40:19] Scissor start...
# ℹ [2026/02/09 15:40:19] Start from raw data...
# ℹ Using "SCT_snn" graph for network.
# ℹ [2026/02/09 15:40:24] Normalizing quantiles of data
# ℹ [2026/02/09 15:41:21] Subsetting data
# ℹ [2026/02/09 15:41:27] Calculating correlation
# -------------------------------------------------------------
# Five-number summary of correlations:
# 0.220832 0.388996 0.414044 0.433624 0.47267
# -------------------------------------------------------------
# ℹ [2026/02/09 15:41:39] Perform cox regression on the given clinical outcomes...
# ✔ [2026/02/09 15:42:43] Statistics data saved to Scissor_inputs.RData.
# ℹ [2026/02/09 15:42:47] Screening...

# ── At alpha = 0.05 ──

# Scissor identified 1789 Scissor+ cells and 1617 Scissor- cells.
# The percentage of selected cell is: 14.533%
# ℹ [2026/02/09 15:45:05] Scissor Ended.
```

### data

The returned structure is a `list`, where the first
slot—`scRNA_data`—contains the Seurat object. Additional slots in the
list store intermediate data generated during the process.

A new column named `Scissor` will be added to the `meta.data` of the
Seurat object, with three possible labels:

- **Positive** denotes cells whose abundance or activity is **positively
  correlated** with the phenotype of interest—specifically, those
  associated with **poor prognosis**.  
- **Negative** denotes cells that are **negatively correlated** with the
  phenotype (i.e., potentially protective or associated with better
  outcomes).  
- **Neutral** cells can be interpreted as **background** or
  **non-informative** cells—those showing little to no association with
  the phenotype.

### files

A new file named `Scissor_inputs.RData` will be created, which contains
the input data for the Scissor algorithm. You can use the intermediate
data for repeated runs to save time when tuning parameters, avoiding the
need to re-run the entire pipeline from scratch. This is an inherent
feature of the `Scissor`.

### note

Other methods (e.g., **scAB**, **scPAS**, etc.) can also be used. In
this vignette, we only introduce and apply **Scissor**. When using
alternative algorithms, remember to explicitly specify `assay = "SCT"`,
as these functions typically default to `assay = "RNA"`, whereas our
data has been processed with `SCTransform`.

## Visualization of Screened Cells

Finally we can visualize the results. Here, we just provide a brief
demonstration using Seurat’s built-in visualization functions, as
visualization preferences vary from user to user. If you need additional
features, feel free to request them in an issue. :)

Let’s first see the spatial position of the Positive cells.

``` r

positive_cell <- colnames(res$scRNA_data)[
  res$scRNA_data$scissor == "Positive"
]

p <- Seurat::SpatialDimPlot(
  res$scRNA_data,
  cells.highlight = positive_cell,
  ncol = 3L,
  image.alpha = 0.5,
  cols.highlight = c("#ff3333", "#CECECE"),
  alpha = c(0.8, 1)
) &
  ggplot2::theme(legend.position = "none")

# ggplot2::ggsave(
#   "vignettes/example_figures/spatial_dim_plot.png",
#   p,
#   width = 14,
#   height = 7
# )
```

[![spatial_dim_plot](example_figures/spatial_dim_plot.png)](NA)

To visualize the cellular composition of Positive cells:

``` r

p2 <- Seurat::SpatialDimPlot(
  spatials,
  group.by = "mllmcelltype_cell_type",
  ncol = 3
)

# ggplot2::ggsave(
#   "vignettes/example_figures/mllmcelltype_spatial.png",
#   p2,
#   width = 14,
#   height = 7
# )
```

[![mllmcelltype_spatial](example_figures/mllmcelltype_spatial.png)](NA)

To statistically assess the composition of Positive cells, you can do
the following:

``` r

table(res$scRNA_data$mllmcelltype_cell_type[
  res$scRNA_data$scissor == "Positive"
])
```

``` r

#     Acinar cells           Adipocytes              B cells         Chondrocytes    Endothelial cells
#              106                  100                   98                  566                   27
#      Enterocytes     Epithelial cells          Fibroblasts        Gastric cells         Goblet cells
#              106                  406                  319                   23                    9
#    Keratinocytes          Macrophages    Mesothelial cells Neuroendocrine cells  Smooth muscle cells
#                7                   17                    2                   43                   54
```

The proportion of Positive cells within each cell type can be calculated
using the following code:

``` r

cell_total <- table(res$scRNA_data$mllmcelltype_cell_type)
cell_positive <- table(res$scRNA_data$mllmcelltype_cell_type[
  res$scRNA_data$scissor == "Positive"
])

ratios <- cell_positive[match(names(cell_total), names(cell_positive))] /
  cell_total
ratios[is.na(ratios)] <- 0

result <- sort(round(ratios, 4), decreasing = TRUE)

result
```

``` r

#     Chondrocytes  Smooth muscle cells         Acinar cells           Adipocytes          Fibroblasts
#           0.6152               0.2126               0.1785               0.1580               0.0940
# Epithelial cells          Enterocytes              B cells        Gastric cells Neuroendocrine cells
#           0.0846               0.0525               0.0424               0.0420               0.0246
#      Macrophages         Goblet cells    Endothelial cells    Mesothelial cells        Keratinocytes
#           0.0193               0.0139               0.0127               0.0083               0.0062
#          T cells
#           0.0000
```

## Sessioninfo

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] Seurat_5.5.0       SeuratObject_5.4.0 sp_2.2-1           SigBridgeR_3.7.0  
#> 
#> loaded via a namespace (and not attached):
#>   [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3         
#>   [4] rlang_1.2.0            magrittr_2.0.5         RcppAnnoy_0.0.23      
#>   [7] otel_0.2.0             spatstat.geom_3.7-3    matrixStats_1.5.0     
#>  [10] ggridges_0.5.7         compiler_4.6.0         png_0.1-9             
#>  [13] systemfonts_1.3.2      vctrs_0.7.3            reshape2_1.4.5        
#>  [16] stringr_1.6.0          pkgconfig_2.0.3        fastmap_1.2.0         
#>  [19] promises_1.5.0         rmarkdown_2.31         ragg_1.5.2            
#>  [22] purrr_1.2.2            xfun_0.57              cachem_1.1.0          
#>  [25] jsonlite_2.0.0         goftest_1.2-3          later_1.4.8           
#>  [28] spatstat.utils_3.2-3   irlba_2.3.7            parallel_4.6.0        
#>  [31] cluster_2.1.8.2        R6_2.6.1               ica_1.0-3             
#>  [34] spatstat.data_3.1-9    bslib_0.11.0           stringi_1.8.7         
#>  [37] RColorBrewer_1.1-3     reticulate_1.46.0      spatstat.univar_3.2-0 
#>  [40] parallelly_1.47.0      lmtest_0.9-40          jquerylib_0.1.4       
#>  [43] scattermore_1.2        Rcpp_1.1.1-1.1         knitr_1.51            
#>  [46] tensor_1.5.1           future.apply_1.20.2    zoo_1.8-15            
#>  [49] sctransform_0.4.3      httpuv_1.6.17          Matrix_1.7-5          
#>  [52] splines_4.6.0          igraph_2.3.1           tidyselect_1.2.1      
#>  [55] abind_1.4-8            yaml_2.3.12            spatstat.random_3.4-5 
#>  [58] spatstat.explore_3.8-0 codetools_0.2-20       miniUI_0.1.2          
#>  [61] listenv_0.10.1         lattice_0.22-9         tibble_3.3.1          
#>  [64] plyr_1.8.9             shiny_1.13.0           S7_0.2.2              
#>  [67] ROCR_1.0-12            evaluate_1.0.5         Rtsne_0.17            
#>  [70] future_1.70.0          fastDummies_1.7.6      desc_1.4.3            
#>  [73] survival_3.8-6         polyclip_1.10-7        fitdistrplus_1.2-6    
#>  [76] pillar_1.11.1          KernSmooth_2.23-26     plotly_4.12.0         
#>  [79] generics_0.1.4         RcppHNSW_0.6.0         ggplot2_4.0.3         
#>  [82] scales_1.4.0           globals_0.19.1         xtable_1.8-8          
#>  [85] glue_1.8.1             lazyeval_0.2.3         tools_4.6.0           
#>  [88] data.table_1.18.4      RSpectra_0.16-2        RANN_2.6.2            
#>  [91] fs_2.1.0               dotCall64_1.2          cowplot_1.2.0         
#>  [94] grid_4.6.0             tidyr_1.3.2            nlme_3.1-169          
#>  [97] patchwork_1.3.2        cli_3.6.6              spatstat.sparse_3.1-0 
#> [100] textshaping_1.0.5      spam_2.11-3            viridisLite_0.4.3     
#> [103] dplyr_1.2.1            uwot_0.2.4             gtable_0.3.6          
#> [106] sass_0.4.10            digest_0.6.39          progressr_0.19.0      
#> [109] ggrepel_0.9.8          htmlwidgets_1.6.4      farver_2.1.2          
#> [112] htmltools_0.5.9        pkgdown_2.2.0          lifecycle_1.0.5       
#> [115] httr_1.4.8             mime_0.13              MASS_7.3-65
```
