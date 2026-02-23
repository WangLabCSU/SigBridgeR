# **SigBridgeR** [![sigbridger website](reference/figures/logo_white.png)](https://wanglabcsu.github.io/SigBridgeR/)

[![Repo
Status](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License:
GPL3](https://img.shields.io/badge/license-GPL3-blue.svg)](https://cran.r-project.org/web/licenses/GPL3)
[![](https://img.shields.io/badge/devel%20version-3.4.0-blue.svg)](https://github.com/WangLabCSU/SigBridgeR)
[![SigBridgeR status
badge](https://wanglabcsu.r-universe.dev/SigBridgeR/badges/version)](https://wanglabcsu.r-universe.dev/SigBridgeR)
[![R CMD
check](https://github.com/WangLabCSU/SigBridgeR/workflows/R-CMD-check/badge.svg)](https://github.com/WangLabCSU/SigBridgeR/actions)
[![registry status
badge](https://wanglabcsu.r-universe.dev/badges/:registry)](https://wanglabcsu.r-universe.dev/)

------------------------------------------------------------------------

## 🌐 Overview

SigBridgeR integrates multiple algorithms, using single-cell RNA
sequencing data, bulk expression data, and sample-related phenotypic
data, to identify the cells most closely associated with the phenotypic
data, performing as a bridge to existing tools.

## 🔧 Installation

Usually we recommend installing the latest release from GitHub because
of the latest features and bug fixes.

1.  Install the development version from GitHub:

``` R
if (!requireNamespace("pak")) {
  install.packages(
    "pak",
    repos = sprintf(
      "https://r-lib.github.io/p/pak/stable/%s/%s/%s",
      .Platform$pkgType,
      R.Version()$os,
      R.Version()$arch
    )
  )
}
pak::pkg_install("WangLabCSU/SigBridgeR")
```

1.  Install from r-universe:

``` R
install.packages("SigBridgeR", repos = "https://wanglabcsu.r-universe.dev")
```

**It is recommended to install the following packages:**

For better performance:

``` R
pak::pkg_install(c(
  # faster computation
  "sparseMatrixStats",
  "matrixStats",
  "preprocessCore",
  "WGCNA",
  "tidyr",
  "matrixTests",
  "KernSmooth",
  "cheapr",
  # better gene symbol conversion
  "scCustomize",
  # parallel computation
  "furrr",
  "future"
))
```

For seamless integration with other file types such as `.h5ad`:

``` R
pak::pkg_install("anndata")
# or
pak::pkg_install("anndataR") # both are supported
```

For visualization:

``` R
pak::pkg_install(c(
  "ggplot2",
  "randomcoloR", # or RColorBrewer
  "ggupset", # for upset plot
  "patchwork", # for fraction stack plot
  "ggforce", # for pca plot
  "ggVennDiagram" # for venn diagram
))
```

To use the built-in cell annotation methods:

``` R
pak::pkg_install(c(
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

To add custom extension functions to SigBridgeR:

``` R
pak::pkg_install(c(
  "tictoc",
  "codetools",
  "knitr",
  "lintr",
  "rstudioapi",
  "yonicd/tidycheckUsage"
))
```

To reproduce the tutorial to learn more usage:

``` R
pak::pkg_install(c(
  "zeallot",
  "here",
  "org.Hs.eg.db",
  "processx"
))
```

## 📓 Documentation

Get Started:

- View [Github Webpage](https://wanglabcsu.github.io/SigBridgeR/)
- [A Quick Started
  Guide](https://wanglabcsu.github.io/sigbridger/vignettes/Quick_Start.md)
- [Full
  Tutorial](https://wanglabcsu.github.io/sigbridger/vignettes/Full_Tutorial.md)
  for more details
- Use `?SigBridgeR::function_name` to access the help documents in R.

If you encounter problems, please check:

- the [Troubleshooting
  Guide](https://wanglabcsu.github.io/sigbridger/vignettes/Troubleshooting.md),
  or
- the [Github issues](https://github.com/WangLabCSU/SigBridgeR/issues)
  page if you want to file bug reports or feature requests

Let us know if you have ideas to make this project better!
