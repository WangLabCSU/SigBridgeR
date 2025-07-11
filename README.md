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

You can install the **SigBridgeR** using the following options:

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
-   [Full Tutorial](vignettes/Full_Tutorial.md) for informed details

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
