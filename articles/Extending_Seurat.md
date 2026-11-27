# Extending SigBridgeR: A Guide for Extend Single Cell Preprocess Pipeline

[![Lifecycle:
deprecated](https://img.shields.io/badge/lifecycle-deprecated-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#deprecated)

> **Deprecated.** This vignette is deprecated as of SigBridgeR 4.0.0 and
> will be removed in a future release.
> [`SCPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcess.md)
> and the registry-based preprocessing pipeline are superseded

``` r

library(SigBridgeR)
#> ✔ SigBridgeR v4.0.0 loaded [471 ms]
```

## Extend Seurat Methods

In general, any Seurat function can be registered for use within
`SCPreProcess`. However, due to the character-length limitation of the
pipeline code string, at most **52** functions can be registered-one for
each uppercase and lowercase letter (A-Z, a-z).

By default, the registry is defined as follows:

**Pipeline Code Table:**

| Code  | Function               | Description                           |
|-------|------------------------|---------------------------------------|
| **o** | `CreateSeuratObject`   | **Required.** Must be the first step. |
| **n** | `NormalizeData`        | Standard normalization.               |
| **s** | `ScaleData`            | Scales data for PCA.                  |
| **v** | `FindVariableFeatures` | Selects highly variable genes.        |
| **p** | `RunPCA`               | Principal Component Analysis.         |
| **e** | `FindNeighbors`        | Computes SNN graph.                   |
| **c** | `FindClusters`         | Louvain algorithm clustering.         |
| **t** | `RunTSNE`              | t-SNE reduction.                      |
| **u** | `RunUMAP`              | UMAP reduction.                       |
| **r** | `SCTransform`          | **SCT workflow.** Replaces n, s, v.   |

For example, if we need to use
[`Seurat::LoadXenium`](https://satijalab.org/seurat/reference/ReadXenium.html)
to load Xenium data for screening, register this function first.

``` r

Register(x = Seurat::LoadXenium)
#> ✔ Registered x to <SCPreProcessStrategy>
```

To check:

``` r

names(SCPreProcessStrategy)
#>  [1] "t" "u" "v" "x" "c" "e" "i" "n" "o" "p" "r" "s"
```

``` r

InterceptStrategy("SCPreProcessStrategy")
#> # A tibble: 12 × 2
#>    name  args      
#>    <chr> <list>    
#>  1 t     <chr [2]> 
#>  2 u     <chr [2]> 
#>  3 v     <chr [2]> 
#>  4 x     <chr [10]>
#>  5 c     <chr [2]> 
#>  6 e     <chr [2]> 
#>  7 i     <chr [8]> 
#>  8 n     <chr [2]> 
#>  9 o     <chr [7]> 
#> 10 p     <chr [1]> 
#> 11 r     <chr [2]> 
#> 12 s     <chr [2]>
```

overwrite:

``` r

Register(x = Seurat::LoadXenium, overwrite = TRUE)
#> ✔ Registered x to <SCPreProcessStrategy>
```

### Session Info

    #> R version 4.6.2 (2026-10-29)
    #> Platform: x86_64-pc-linux-gnu
    #> Running under: Ubuntu 26.04.1 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
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
    #> [1] SigBridgeR_4.0.0
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] tidyselect_1.2.1      dplyr_1.2.1           farver_2.1.2         
    #>  [4] S7_1.0.0              fastmap_1.2.0         RANN_2.6.3           
    #>  [7] digest_0.6.39         dotCall64_1.2         lifecycle_1.0.5      
    #> [10] cluster_2.1.8.3       SeuratObject_5.4.0    magrittr_2.0.5       
    #> [13] compiler_4.6.2        rlang_1.3.0           sass_0.5.0           
    #> [16] tools_4.6.2           utf8_1.2.6            igraph_2.3.4         
    #> [19] yaml_2.3.12           data.table_1.18.6.1   knitr_1.52           
    #> [22] htmlwidgets_1.6.4     sp_2.2-4              reticulate_1.47.0    
    #> [25] plyr_1.8.9            RColorBrewer_1.1-3    KernSmooth_2.23-27   
    #> [28] Rtsne_0.17            purrr_1.2.2           desc_1.4.3           
    #> [31] grid_4.6.2            RcppHNSW_0.7.0        fastDummies_1.7.6    
    #> [34] future_1.82.0         progressr_1.1.0       ggplot2_4.0.4        
    #> [37] globals_0.20.0        scales_1.4.0          MASS_7.3-66          
    #> [40] ggridges_0.5.7        cli_3.6.6             rmarkdown_2.32       
    #> [43] ragg_1.5.2            generics_0.1.4        Seurat_6.0.1         
    #> [46] otel_0.3.0            ica_1.0-3             future.apply_1.21.0  
    #> [49] RSpectra_0.16-2       sctransform_0.4.3     httr_1.4.9           
    #> [52] reshape2_1.4.5        pbapply_1.7-5         cachem_1.1.0         
    #> [55] stringr_1.6.0         parallel_4.6.2        matrixStats_1.5.0    
    #> [58] vctrs_0.7.3           Matrix_1.7-6          jsonlite_2.0.0       
    #> [61] patchwork_1.3.2       ggrepel_0.9.8         irlba_2.4.1          
    #> [64] scattermore_1.2       listenv_1.1.0         systemfonts_1.3.4    
    #> [67] jquerylib_0.1.4       pls_2.9-0             glue_1.8.1           
    #> [70] parallelly_1.50.0     spam_2.11-4           pkgdown_2.2.1        
    #> [73] chk_0.11.0            SigBridgeRUtils_0.2.7 codetools_0.2-20     
    #> [76] uwot_0.2.5            cowplot_1.2.0         stringi_1.8.9        
    #> [79] RcppAnnoy_0.0.23      gtable_0.3.6          tibble_3.3.1         
    #> [82] pillar_1.11.1         htmltools_0.5.9       R6_2.6.1             
    #> [85] textshaping_1.1.0     evaluate_1.0.5        lattice_0.23-1       
    #> [88] png_0.1-9             bslib_0.12.0          Rcpp_1.1.2           
    #> [91] gridExtra_2.3.1       xfun_0.62             fs_2.1.0             
    #> [94] pkgconfig_2.0.3
