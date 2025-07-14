# **SigBridgeR** Under development

[![CRAN Status](https://www.r-pkg.org/badges/version/SigBridgeR)](https://cran.r-project.org/package=SigBridgeR) [![R-CMD-check](https://github.com/WangLabCSU/SigBridgeR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/WangLabCSU/SigBridgeR/actions) [![Codecov test coverage](https://codecov.io/gh/WangLabCSU/SigBridgeR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/WangLabCSU/SigBridgeR) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

An R package for the integration of single-cell RNA-seq data, mutational signatures and bulk expression data.

------------------------------------------------------------------------

## 🌐 Overview

SigBridgeR integrates multiple algorithms to identify the cells most closely associated with the Mutational Signatures phenotype, preforming as an birdge to the existing tools.

SigbridgeR integrates algorithms from the following packages:

-   [Github-sunduanchen/Scissor](https://github.com/sunduanchen/Scissor)
-   [Github-Qinran-Zhang/scAB](https://github.com/Qinran-Zhang/scAB/)
-   [Github-WangX-Lab/ScPP](https://github.com/WangX-Lab/ScPP)
-   [Github-aiminXie/scPAS](https://github.com/aiminXie/scPAS)

## 🔧 Installation

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
remotes::install_github("WangLabCSU/SigBridgeR")
```

## 📓 Documentation

Use `?SigBridgeR` to access the help documents in  R.

Key resources:

-   [Quick Started Guide](vignettes/Quick_Start.md) 
-   [Full Tutorial](vignettes/Full_Tutorial.md) for more details

If you encounter problems, please see:

-   [Troubleshooting Guide](vignettes/Troubleshooting.md)

Other information:

- What is *Mutational Signature*? 
  - [COSMIC | Mutational Signatures - An Introduction ](https://cancer.sanger.ac.uk/cosmic)
  - [View in Wiki](https://en.wikipedia.org/wiki/Mutational_signatures)
  - [Nature | Mutational signatures: emerging concepts, caveats and clinical applications](https://www.nature.com/articles/s41568-021-00377-7)
- What is *Single Cell Sequencing*? 
  - [Veiw in Wiki](https://en.wikipedia.org/wiki/Single-cell_sequencing)
- What is *RNA-seq*? 
  - [View in Wiki](https://en.wikipedia.org/wiki/RNA-Seq)

## 📮 Contact

For support or questions:

Maintainer: Exceret
