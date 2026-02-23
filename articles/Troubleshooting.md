# Troubleshooting

### Overview

This vignette provides solutions to common problems you might encounter
when using this package. If you cannot find your problem here, please
file an issue on our GitHub repository.

  

------------------------------------------------------------------------

## Scissor

> **Error in normalize.quantiles(dataset0) : ERROR; return code from
> pthread_create() is 22**

To solve this problem, you need to install `preprocessCore` without
threading support. Try:

``` shell
# ! shell
git clone https://github.com/bmbolstad/preprocessCore.git
cd preprocessCore
R CMD INSTALL --configure-args="--disable-threading"  .
```

or

``` r
BiocManager::install(
  "preprocessCore",
  configure.args = "--disable-threading",
  force = TRUE
)
```

See
[bioconductor_docker/issues/22](https://github.com/Bioconductor/bioconductor_docker/issues/22),
[Scissor/issues/15](https://github.com/sunduanchen/Scissor/issues/15)
for more details.

  

------------------------------------------------------------------------

> **Error at alpha=0.05:subscript out of bounds**
>
> **Error in `Scissor.v5.optimized()`:**
>
> **! object ‘fit0’ not found**

This may be due to two reasons: first, a mismatch in the dimensions of
the bulk expression data and the phenotype data; second, incorrect
column names in the survival phenotype data leading to a failure to
match.

To check dimension mismatch：

``` r
ncol(your_bulk_data) == nrow(your_phenotype_data) # should be TRUE

all(unique(colnames(your_bulk_data)) == unique(rownames(your_bulk_data))) # should be TRUE

all(order(colnames(your_bulk_data)) == order(rownames(your_phenotype_data))) # should be TRUE
```

Survival phenotype column names should be formatted as `time` and
`status`, ensuring correct capitalization and spelling:

``` r
head(survival_phenotype) # case-sensitive
#           time status
# GSM70130 34.80      0
# GSM70131 35.67      0
# GSM70136 43.37      0
# GSM70138 60.77      0
# GSM70140 33.80      1
# GSM70144 58.53      0
```

  

------------------------------------------------------------------------

## Package

> **Error in `function_name()`:**
>
> **! lazy-load database
> ‘/home/user/R/x86_64-pc-linux-gnu-library/4.4/SigBridgeR/R/SigBridgeR.
> rdb’ is corrupt**

An error occurred during installation, causing the package to be
corrupted. Try reinstalling the package:

``` r
# I reccomend restarting R/RStudio before reinstalling
remove.packages("SigBridgeR")

detach("package:SigBridgeR", unload = TRUE)

if (!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("WangLabCSU/SigBridgeR")
```

  

------------------------------------------------------------------------

## Visualization

> **Error:**
>
> **! Invalid syntax: ‘c(scissor_umap, scpas_umap)’**

This may be due to the R environment being contaminated, which prevents
the use of `%<-%`. You can try restarting the R environment or clean up
the environment.

``` r
# restart R/RStudio
.rs.api.restartSession()
```

``` r
rm(list = ls(all.names = TRUE))
```

  

------------------------------------------------------------------------

#### scPP

> **Error:**
>
> **! Detected n gene(s) with zero variance:**
>
> **ℹ “gene name(s)”**

This is due to the presence of genes with zero variance in the bulk
expression data when you are using `scPP` and `binary` phenotype. This
indicates that the expression levels of one (or several) genes are
nearly identical across different samples. You should check your data.

  

------------------------------------------------------------------------

> **✖ $$2025/09/2308:53:51$$ Fewer than 20% of the genes in the gene
> sets are included in the rankings.Check wether the gene IDs in the
> ‘rankings’ and ‘geneSets’ match.**
>
> **ℹ $$2025/09/2308:53:51$$ scPP screening exit.**

This issue arises from the single-cell processing, which filtered out
too many genes and cells. Consider Adjusting `min.cells` and
`min.features` to a smaller value.:

  

#### DEGAS

------------------------------------------------------------------------

> **Warning in info\$envs : partial match of ‘envs’ to ‘envs
> directories’**
>
> **Error in
> [`reticulate::use_condaenv()`](https://rstudio.github.io/reticulate/reference/use_python.html):**
>
> **! Unable to locate conda environment ‘r-reticulate-degas’.**

`reticulate` cannot find the relevant environment because a Python
backend has been registered. One solution is to pass in the Python path
of the environment instead of the environment name.

``` r
envs <- ListPyEnv()
head(envs) # goes like this
#                                                            name                                                  python  type
# /home/user/miniconda3                                       base                         /home/user/miniconda3/bin/python conda
# /home/user/miniconda3/envs/r-reticulate-degas r-reticulate-degas /home/user/miniconda3/envs/r-reticulate-degas/bin/python conda
```

For example, if wanna use the `r-reticulate-degas` environment, pass its
Python location to
[`reticulate::use_condaenv()`](https://rstudio.github.io/reticulate/reference/use_python.html).

``` r
py_path <- envs[envs$name == "r-reticulate-degas", "python"]
# /home/user/miniconda3/envs/r-reticulate-degas/bin/python
reticulate::use_condaenv(py_path)
```

  

------------------------------------------------------------------------

> **Traceback (most recent call last):**
>
> **File
> “/home/user/R/Project/R_code/SigBridgeR/Tmp/tmp/BlankClassMTL.py”,
> line 1, in**
>
> **import tensorflow as tf \#NEED**
>
> **ModuleNotFoundError: No module named ‘tensorflow’**
>
> **Warning in file(file, “r”) :**
>
> **cannot open file ‘tmp//Activations.csv’: No such file or directory**

This is due to the issue with the tensorflow package. Use the following
method to check if the tensorflow package can be imported.

``` python
# make sure you are using the correct python, i.e. r-reticulate-degas
import sys

sys.executable  # in an environment, it should be something like /home/user/miniconda3/envs/r-reticulate-degas/bin/python3

import tensorflow as tf

tf.__version__
```

``` shell
# make sure you've entered the Conda environment you want to check, i.e. r-reticulate-degas
conda list tensorflow

# packages in environment at /home/user/miniconda3/envs/r-reticulate-degas:                                                               
#                                                                                                                                        
# Name                     Version          Build               Channel                                                                  
# tensorflow                 2.4.1            mkl_py39h4683426_0
# tensorflow-base            2.4.1            mkl_py39h43e0292_0
# tensorflow-estimator       2.6.0            py39he80948d_0      https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
```

If the output is `ModuleNotFoundError: No module named 'tensorflow'`,
you need to install the tensorflow package, and strictly control the
version of the package (version 2.4.1 has been tested and is feasible)

  

------------------------------------------------------------------------

#### Seurat

> **Excessive time during running SCTransfomm**

Output message is probably like this:

``` r
# Running SCTransform on assay: Xenium
# vst.flavor='v2' set. Using model with fixed slope and excluding poisson genes.
# Calculating cell attributes from input UMI matrix: log_umi
# Variance stabilizing transformation of count matrix of size 248 by 36553
# Model formula is y ~ log_umi
# Get Negative Binomial regression parameters per gene
# Using 248 genes, 5000 cells
# Second step: Get residuals using fitted parameters for 248 genes
# Computing corrected count matrix for 248 genes
# Calculating gene attributes
# Wall clock passed: Time difference of 7.166213 secs
# Determine variable features
# ! stuck here.
```

This is a bug from Seurat, see [Seurat
\#10153](https://github.com/satijalab/seurat/issues/10153)

Try update `Seurat` to latest version

## Other

> **Why does SigBridgeR hardcode reticulate::use_python() internally?
> This contradicts the recommended practices of the R community**

We are indeed aware that this practice deviates from the expectations of
the R community. However, given that multiple algorithms embedded within
the package require distinct Python installations and environments—and
considering that `reticulate` binds a single Python environment to the
entire R session without allowing runtime switching—each screening
algorithm would otherwise necessitate a separate R session to operate
correctly.

Meanwhile, wrapping functions in subprocesses would increase usability
and learning overhead for users, while invoking Python interpreters via
system commands would require additional read/write operations to disk
for R-side data, which is inefficient for large-scale data processing.

After careful consideration, we encapsulated
[`reticulate::use_python()`](https://rstudio.github.io/reticulate/reference/use_python.html)
within isolated functions to emulate subprocess-like behavior—allowing
localized Python environment specification—while avoiding interference
with the globally registered Python configuration in the R session.
