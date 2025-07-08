# **SigBridgeR**

[![CRAN Status](https://www.r-pkg.org/badges/version/SigBridgeR)](https://cran.r-project.org/package=SigBridgeR) [![R-CMD-check](https://github.com/yourusername/SigBridgeR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yourusername/SigBridgeR/actions) [![Codecov test coverage](https://codecov.io/gh/yourusername/SigBridgeR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/yourusername/SigBridgeR) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

An R package for the integration of single-cell RNA-seq data, mutational signatures and bulk expression data.

------------------------------------------------------------------------

## Overview

SigBridgeR integrates multiple algorithms to identify the cells most closely associated with the Mutational Signature phenotype, thereby analyzing their roles and functions within the tumor immune microenvironment.

SigbridgeR integrates algorithms from the following packages:

-   [sunduanchen/Scissor](https://github.com/sunduanchen/Scissor)
-   [Qinran-Zhang/scAB](https://github.com/Qinran-Zhang/scAB/)
-   [WangX-Lab/ScPP](https://github.com/WangX-Lab/ScPP)
-   [aiminXie/scPAS](https://github.com/aiminXie/scPAS)

## Installation

### Dependencies

You can use this function to quickly check if all dependencies are installed and if their versions are correct:

```{r check_dependencies}
CheckPkgs <- function(packages) {
  CheckSinglePkg <- function(pkg_spec) {
    installed <- requireNamespace(pkg_spec$pkg, quietly = TRUE)
    current_version <- NA_character_
    
    if (installed) {
      current_version <- as.character(utils::packageVersion(pkg_spec$pkg))
      
      if (!is.null(pkg_spec$version)) {
        installed <- current_version >= pkg_spec$version
      }
    }
    
    data.frame(
      Package = pkg_spec$pkg,
      Required_Version = if (is.null(pkg_spec$version)) "Any" else pkg_spec$version,
      Installed = installed,
      Current_Version = if (installed) current_version else NA,
      stringsAsFactors = FALSE
    )
  }
  
  result <- do.call(rbind, lapply(packages, CheckSinglePkg))
  rownames(result) <- NULL  
    
  return(result)
}

CheckPkgs(list(
  list(pkg = "Seurat", version = "5.0.0"),
  list(pkg = "dplyr"),
  list(pkg = "cli"),
  list(pkg = "AUCell"),
  list(pkg = "future"),
  list(pkg = "ggrepel"),
  list(pkg = "grid"),
  list(pkg = "IDConverter"),
  list(pkg = "Matrix"),
  list(pkg = "patchwork"),
  list(pkg = "scAB"),
  list(pkg = "Scissor"),
  list(pkg = "scPAS"),
  list(pkg = "ScPP"),
  list(pkg = "tibble"),
  list(pkg = "tidyr"))
)

```

Then you can install the **SigBridgeR** using the following options:

### Stable release from CRAN

```{r install_from_cran}
install.packages("SigBridgeR")
```

### Development version from GitHub

```{r install_from_github}
if(!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("yourusername/SigBridgeR")
```

## Documentation

In the R terminal, please use the command `?SigBridgeR` to access the help documents.

Key resources:

-   [Quick Started Guide](vignettes/Quick_Start.md)
-   [Full Tutorial](vignettes/Full_Tutorial.md)

If you encounter problems, please see:

-   [Troubleshooting Guide](vignettes/Troubleshooting.md)

## Function Overview

### Data Processing

-   `SCPreProcess()`:
-   `MSPreProcess()`:
-   `BulkPreProcess()`:
-   `MatchSample()`:

### Core Methods

-   `Screen()`:

### Output Utilities

-   `FetchUMAP()`:

## Contact

For support or questions:\
Maintainer: Exceret\
