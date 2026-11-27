# Extending SigBridgeR: A Guide for Extend SCAnnotate

``` r

library(SigBridgeR)
#> ✔ SigBridgeR v4.0.0 loaded [477 ms]
```

## Extend Cell Type Annotation Methods

By default, SigBridgeR includes three built-in cell type annotation
algorithms: **SingleR**, **CellTypist**, and **mLLMCelltype**-all of
which require additional installation and dependencies (See
[README](https://wanglabcsu.github.io/SigBridgeR/index.html)). They can
all be invoked through the unified interface function `SCAnnotate`.

Extending cell annotation method is easy:

``` r

my_method <- function(seurat_obj, ...) {
  # TODO: Annotate cell type
  return(seurat_obj)
}
my_sc_anno <- AnnotationMethod(
  method_name = "my_method",
  executor = my_method
)
Register(
  my_method = my_sc_anno
)
#> ✔ Registered my_method to <SCAnnotateStrategy>
```

To check:

``` r

names(SCAnnotateStrategy)
#> [1] "mLLMCelltype" "CellTypist"   "SingleR"      "my_method"
```

``` r

InterceptStrategy("SCAnnotateStrategy")
#> # A tibble: 4 × 5
#>   sigbridger_version method_name  method_version executor  name        
#>   <chr>              <chr>        <chr>          <list>    <chr>       
#> 1 4.0.0              mLLMCelltype 2.0.8          <chr [2]> mLLMCelltype
#> 2 4.0.0              CellTypist   1.7.1          <chr [2]> CellTypist  
#> 3 4.0.0              SingleR      2.14.2         <chr [2]> SingleR     
#> 4 4.0.0              my_method    NA             <chr [2]> my_method
```

Here are some recommendations for user-defined cell annotation
functions.

1.  The first input and the first output should both be Seurat objects.
2.  If multiple objects need to be returned, it is recommended to use a
    `list` format. And the first element of the list should be the
    Seurat object.

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
    #>  [1] Matrix_1.7-6          gtable_0.3.6          jsonlite_2.0.0       
    #>  [4] dplyr_1.2.1           compiler_4.6.2        tidyselect_1.2.1     
    #>  [7] Rcpp_1.1.2            jquerylib_0.1.4       png_0.1-9            
    #> [10] systemfonts_1.3.4     scales_1.4.0          textshaping_1.1.0    
    #> [13] yaml_2.3.12           fastmap_1.2.0         lattice_0.23-1       
    #> [16] reticulate_1.47.0     ggplot2_4.0.4         R6_2.6.1             
    #> [19] SigBridgeRUtils_0.2.7 generics_0.1.4        knitr_1.52           
    #> [22] htmlwidgets_1.6.4     tibble_3.3.1          desc_1.4.3           
    #> [25] bslib_0.12.0          pillar_1.11.1         RColorBrewer_1.1-3   
    #> [28] rlang_1.3.0           utf8_1.2.6            cachem_1.1.0         
    #> [31] xfun_0.62             fs_2.1.0              sass_0.5.0           
    #> [34] S7_1.0.0              otel_0.3.0            cli_3.6.6            
    #> [37] pkgdown_2.2.1         magrittr_2.0.5        digest_0.6.39        
    #> [40] grid_4.6.2            lifecycle_1.0.5       vctrs_0.7.3          
    #> [43] evaluate_1.0.5        glue_1.8.1            data.table_1.18.6.1  
    #> [46] farver_2.1.2          ragg_1.5.2            rmarkdown_2.32       
    #> [49] purrr_1.2.2           pkgconfig_2.0.3       tools_4.6.2          
    #> [52] htmltools_0.5.9
