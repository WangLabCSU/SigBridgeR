<!-- README.md is generated from README.Rmd. Please edit that file -->

# SigBridgeR <a href="https://wanglabcsu.github.io/SigBridgeR/"><img src="man/figures/logo_white.png" alt="sigbridger website" align="right" height="139"/></a>

<!-- badges: start -->

[![Project_Status:\_Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) [![Repo_Status](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable) [![License:GPL3](https://img.shields.io/badge/license-GPL3-blue.svg)](https://cran.r-project.org/web/licenses/GPL3) [![Devel_version](https://img.shields.io/badge/devel%20version-4.0.0-blue.svg)](https://github.com/WangLabCSU/SigBridgeR) [![R_CMD_check](https://github.com/WangLabCSU/SigBridgeR/workflows/R-CMD-check/badge.svg)](https://github.com/WangLabCSU/SigBridgeR/actions) [![Ask_DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/WangLabCSU/SigBridgeR) [![Codecov_test_coverage](https://codecov.io/gh/WangLabCSU/SigBridgeR/graph/badge.svg)](https://app.codecov.io/gh/WangLabCSU/SigBridgeR) [![registry_status_badge](https://wanglabcsu.r-universe.dev/badges/:registry)](https://wanglabcsu.r-universe.dev/) [![code_size](https://img.shields.io/github/languages/code-size/https://github.com/WangLabCSU/SigBridgeR.svg)](https://github.com/WangLabCSU/SigBridgeR) <!-- badges: end -->

## 🌐 Overview

SigBridgeR integrates multiple algorithms, using single-cell RNA sequencing data, bulk expression data, and sample-related phenotypic data, to identify the cells most closely associated with the phenotypic data, performing as a bridge to existing tools.

<p align="center">

<img src="SigBridgeR_workflow.png" />

</p>

**Current Status**

|  Algorithm  |          Supported Phenotypes          | Cache Support |  License   |
|:-----------:|:--------------------------------------:|:-------------:|:----------:|
| **Scissor** | `binary`<br>`survival`<br>`continuous` |      ✅       |   GPL-3    |
|  **scPAS**  | `binary`<br>`survival`<br>`continuous` |               |   GPL-3    |
|  **scAB**   |         `binary`<br>`survival`         |      ✅       |   GPL-3    |
|  **scPP**   | `binary`<br>`survival`<br>`continuous` |               |   GPL-3    |
|  **DEGAS**  | `binary`<br>`survival`<br>`continuous` |               |    MIT     |
| **LP-SGL**  | `binary`<br>`survival`<br>`continuous` |               | No License |
|  **PIPET**  | `binary`<br>`survival`<br>`continuous` |               |   GPL-3    |
| **SIDISH**  |               `survival`               |               |    MIT     |
| **SCIPAC**  | `binary`<br>`survival`<br>`continuous` |               |    MIT     |
| **TiRank**  | `binary`<br>`survival`<br>`continuous` |               |    MIT     |

## 🔧 Installatation

Usually we recommend installing the latest release from GitHub because of the latest features and bug fixes.

1.  Install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("WangLabCSU/SigBridgeR")
```

2.  Install from r-universe:

``` r
install.packages("SigBridgeR", repos = "https://wanglabcsu.r-universe.dev")
```

**unnecessary but recommended**:

<details>

<summary>
For better performance:
</summary>

``` r
pak::pak(
  # faster computation
  "sparseMatrixStats",
  "matrixStats",
  "preprocessCore",
  "tidyr",
  "matrixTests",
  "KernSmooth",
  "cheapr",
  # better gene symbol conversion
  "scCustomize",
  "WGCNA"
))
```

</details>

<details>

<summary>
For seamless integration with single-cell RNA-seq data stored in `.h5ad`:
</summary>

``` r
pak::pak("anndata")
# or
pak::pak("anndataR") # both are supported
```

</details>

<details>

<summary>
For visualization:
</summary>

``` r
pak::pak(c(
  "ggplot2",
  "randomcoloR", # or RColorBrewer
  "ggupset", # for upset plot
  "patchwork", # for fraction stack plot
  "ggforce", # for pca plot
  "ggVennDiagram" # for venn diagram
))
```

</details>

<details>

<summary>
To use the built-in cell annotation methods:
</summary>

``` r
pak::pak(c(
  # SingleR
  "SingleR-inc/SingleR",
  "celldex",
  # mLLMCelltype
  "mLLMCelltype",
  "plyr",
  # CellTypist
  "reticulate",
  "AnnDataR"
))
```

</details>

<details>

<summary>
To add custom extension functions to SigBridgeR:
</summary>

``` r
pak::pak(c(
  "tictoc",
  "codetools",
  "lintr",
  "rstudioapi"
))
```

</details>

<!-- 
## 🤝 Community
&#10;SigBridgeR allows the community developers to integrate new algorithms without modifying the core package. The following community-developed plugins are maintained by their respective authors.
&#10;
``` r
# Install community plugins via pak:
pak::pak(c(
  # Example: "username/sigbridger-plugin-example"
))
```
&#10;| Name | Author | Status | Repo Link |
|------|--------|--------|-----------|
| *Coming soon* | Be the first to contribute! | - | - |
&#10;To contribute your own plugin, refer to the
[vignette on extending SigBridgeR](https://wanglabcsu.github.io/SigBridgeR/articles/extending.html).
-->

## 📓 Get Started

- View [Github Webpage](https://wanglabcsu.github.io/SigBridgeR/)
- Use `?SigBridgeR::function_name` to access the help documents in R.

If you encounter problems, please check:

- the [Troubleshooting Guide](vignettes/Troubleshooting.md), or
- the [Github issues](https://github.com/WangLabCSU/SigBridgeR/issues) page if you want to file bug reports or feature requests

Let us know if you have ideas to make this project better. Pull requests are welcome!

## 📚 Citation

If you use SigBridgeR in your work, please cite the corresponding papers:

``` r
citation("SigBridgeR")
#> Please cite the following publications when using this package.
#> 
#>   Yang Y, Wang S (2026). _SigBridgeR: Integrative Toolkit for Linking Phenotypes to Cell Subpopulations via scRNA-seq and
#>   Bulk Data_. R package version 4.0.0, <https://wanglabcsu.github.io/sigbridger>.
#> 
#>   Sun et al. Identifying phenotype-associated subpopulations by integrating bulk and single-cell sequencing data. Nat
#>   Biotechnol (2022) [Scissor]
#> 
#>   Johnson et al. Diagnostic Evidence GAuge of Single cells (DEGAS): a flexible deep transfer learning framework for
#>   prioritizing cells in relation to disease. Genome Med (2022) [DEGAS]
#> 
#>   Zhang et al. scAB detects multiresolution cell states with clinical significance by integrating single-cell genomics
#>   and bulk sequencing data. Nucleic Acids Res (2022) [scAB]
#> 
#>   Li et al. Identifying phenotype-associated subpopulations through LP_SGL. Brief Bioinform (2024) [LP_SGL]
#> 
#>   Gan et al. SCIPAC: quantitative estimation of cell-phenotype associations. Genome Biol (2024) [SCIPAC]
#> 
#>   Ruan et al. PIPET: predicting relevant subpopulations in single-cell data using phenotypic information from bulk data.
#>   Brief Bioinform (2024) [PIPET]
#> 
#>   Xie et al. scPAS: single-cell phenotype-associated subpopulation identifier. Brief Bioinform (2025) [scPAS]
#> 
#>   Jolasun et al. SIDISH integrates single-cell and bulk transcriptomics to identify high-risk cells and guide precision
#>   therapeutics through in silico perturbation. Nat Commun (2025) [SIDISH]
#> 
#>   He et al. Inferring Phenotypes of Single Cells Based on the Expression Profiles of Phenotype-Associated Marker Genes in
#>   Bulks and Single Cells. Interdiscip Sci (2026) [ScPP]
#> 
#>   Lin et al. TiRank prioritizes phenotypic niches in tumor microenvironment for clinical biomarker discovery. Genome Med
#>   (2026) [TiRank]
#> 
#> To see these entries in BibTeX format, use 'print(<citation>, bibtex=TRUE)', 'toBibtex(.)', or set
#> 'options(citation.bibtex.max=999)'.
```

## 🗺️ Similar Projects

[scSurvival](https://github.com/cliffren/scSurvival): Single-cell expression matrix (log-normalized + HVG-selected) + survival data (optional clinical covariates and batch labels) -\> Survival-associated cell subpopulations

[CellPhenoX](https://github.com/fanzhanglab/pyCellPhenoX): Single-cell multi-omics data + bulk-level clinical variables, covariates (optional interaction effect terms) -\> interpretable score per cell

[scPrognosis](https://github.com/XiaomeiLi1/scPrognosis): scRNA-seq count matrix (imputed by MAGIC + filtered for low coverage/expression) + bulk RNA-seq expression matrix (with matched survival time and event status) -\> breast cancer prognostic gene signatures and Cox PH risk prediction model

[SCellBOW](https://github.com/cellsemantics/SCellBOW): source scRNA-seq expression matrix + target scRNA-seq expression matrix + (optional) bulk RNA-seq expression matrix with paired survival data -\> cell embeddings, cluster assignments, UMAP visualizations, and phenotype-algebra-derived risk scores with survival probability curves for individual cell subpopulations.

[scPER](https://github.com/BrianLlll/scPER): Single-cell RNA-seq data + bulk RNA-seq data + celltype annotation (optional batch labels) -\> phenotype-associated cell populations
