# Full Tutorial

## Full Tutorial for SigBridgeR

### 0. Preface

SigBridgeR (short for **Sig**nificant cell-to-phenotype **Bridge** in
**R**) is an R package for screening cells highly associated with
phenotype data using single-cell RNA-seq, bulk RNA expression and sample
related phenotype data (e.g. patient survival, age, etc). It integrates
many single cell phenotypic screening methods ([8.
References](#id_8-references)) and provides unified preprocessing,
parameter tuning, and visualization approaches,performing as a unified
integration panel.

------------------------------------------------------------------------

### 1. Load and Preprocess data

First load the packages, you will see a version number message
indicating successful loading:

``` r

library(SigBridgeR)
```

    ## ✔ SigBridgeR v4.0.0 loaded [281 ms]

#### 1.1 Single-cell RNA-seq Data

You need to preprocess your single-cell RNA-seq data into a Seurat
object. We provide some function helping you do this purpose.

##### 1.1.1 Find the Optimal Normalization Method for Seurat

A simple scoring function for selecting a normalization method.

Evaluates 3 dimensions of preprocessing quality:

- **Variance stabilization**: Decoupling of mean-variance relationship
  in normalized expression (lower correlation = better).

- **Biological signal retention**: Preservation of known marker genes
  within highly variable genes (higher retention = better).

- **Dropout robustness**: Removal of technical dropout bias from
  normalized values (lower correlation with dropout rate = better).

Exmaple:

``` r

obj <- Seurat::CreateSeuratObject(counts)

sct <- Seurat::SCTransform(obj)
lognorm <- Seurat::NormalizeData(obj) |>
  Seurat::ScaleData() |>
  Seurat::FindVariableFeatures()

res <- ChooseNormalization(
  sct = sct,
  lognorm = lognorm,
  # * More can be added here
  subset_size = integer(),
  known_hvgs = list(),
  n_hvgs = 2000L,
  low_expressed_thresh = 0.2,
  weight = c(
    variance_stability = 0.4,
    marker_signal = 0.35,
    dropout_robustness = 0.25
  )
)
```

``` r

res$metrics |> format()
res$recommendation
```

Supports simultaneously adding Seurat objects preprocessed with multiple
methods, ensuring they are equivalent to the standard workflow of
`NormalizeData` + `ScaleData` + `FindVariableFeatures`.

**Parameters**

- `subset_size`: Size of the single-cell data subset (i.e., number of
  cells) to be evaludated.
- `known_hvgs`: list or character. Known highly variable genes, provided
  as a reference baseline.
- `n_hvgs`: Maximum number of highly variable genes to consider.
- `low_expressed_thresh`: Percentage of genes excluded from
  consideration, ranked by expression level.
- `weight`: Scoring weight that sum to 1.

##### 1.1.2 Finding the Optimal Number of Principle Components

Usually the number if principle components (PCs) is manually set to 10
or 20 according to the elbow plot. However, it is not always the case
that the PCs with the highest variance are the most informative. Here is
a function to help you find the optimal number of PCs.

``` r

ndims <- FindRobustElbow(
  obj = seurat_obj,
  verbose = TRUE,
  ndims = 50
)
```

``` r

ndims
```

This function will draw an elbow plot with each method and the
recommended number of PCs. You can also use the `ndims` parameter to
specify the maximum number of PCs to be tested. The default value is
`50`.

##### 1.1.3 Add miscellaneous information to the Seurat object

SigBridgeR uses
[`AddMisc()`](https://wanglabcsu.github.io/sigbridger/reference/AddMisc.md)
to record what data features or evidence the various screening
algorithms are based on during execution. As a substitute for the
[`SeuratObject::Misc()`](https://satijalab.github.io/seurat-object/reference/Misc.html)

- [`AddMisc()`](https://wanglabcsu.github.io/sigbridger/reference/AddMisc.md)
  : Add miscellaneous information to the Seurat object. Support for
  adding multiple attributes to the `SeuratObject@misc` slot
  simultaneously.

``` r

# basic usage
seurat_obj <- AddMisc(seurat_obj, "QC_stats" = "^MT-")

# Auto-incrementing example when `cover` set to FALSE
seurat_obj <- AddMisc(seurat_obj, markers = c("CD3D", "CD4", "CD8A"))
seurat_obj <- AddMisc(
  seurat_obj,
  markers = c("CD3D", "CD4", "CD8A"),
  cover = FALSE
)

# Add multiple attributes to the `SeuratObject@misc` slot simultaneously
seurat_obj <- AddMisc(seurat_obj, markers1 = "CD3D", markers2 = "CD4")

# Add multiple attributes (stored as a list element) to the `SeuratObject@misc`
seurat_obj <- AddMisc(seurat_obj, list(attr1 = "value1", attr2 = "value2"))
```

``` r

seurat_obj@misc |> str()
```

##### 1.1.4 Add Gene-level Metadata to the Seurat object

SigBridgeR uses
[`AddMetaFeature()`](https://wanglabcsu.github.io/sigbridger/reference/AddMetaFeature.md)
to update the `SeuratObject@assays$RNA@meta.data`

Add a new column to the metadata

``` r

gene_type <- paste0("test", 1:nrow(seurat_obj))

seurat_obj <- AddMetaFeature(seurat_obj, "gene_type" = gene_type)
```

``` r

seurat_obj@assays$RNA@meta.data |> head()
```

Add a data.frame to the metadata, with the column names as the metadata
names.

``` r

seurat_obj <- AddMetaFeature(
  seurat_obj,
  data.frame(gene_type = gene_type, gene_name = rownames(seurat_obj))
)
```

``` r

seurat_obj@assays$RNA@meta.data |> head()
```

Add to different assays

``` r

seurat_obj <- AddMetaFeature(seurat_obj, "gene_type" = gene_type, assay = "RNA")
seurat_obj <- AddMetaFeature(
  seurat_obj,
  "gene_type" = gene_type,
  assay = "ATAC"
)
```

If duplicate column names are detected, they will be suffixed with an
underscore and a number (e.g., `_1`, `_2`) for disambiguation.

##### 1.1.5 Integration of Single Cell Data

**Matrix Integration**

When passing raw matrices, the function performs a “join” operation
based on the union of all genes. Missing values are filled with NA.

``` r

# Create dummy matrices
mat1 <- matrix(
  rpois(50, 5),
  nrow = 10,
  dimnames = list(paste0("G", 1:10), paste0("C", 1:5))
)
mat2 <- matrix(
  rpois(60, 6),
  nrow = 10,
  dimnames = list(paste0("G", 5:14), paste0("C", 1:6))
)

# Integrate matrices with custom prefixes
integrated_mat <- SCIntegrate(BatchA = mat1, BatchB = mat2)
# default prefixed
integrated_mat2 <- SCIntegrate(mat1, mat2)
```

``` r

integrated_mat[1:6, 1:7]
```

    ## 6 x 7 sparse Matrix of class "dgCMatrix"
    ##    BatchA_C1 BatchA_C2 BatchA_C3 BatchA_C4 BatchA_C5 BatchB_C1 BatchB_C2
    ## G1         2         8         4         1         6         .         .
    ## G2         7         3         6         3         5         .         .
    ## G3         5         1         6         4         6         .         .
    ## G4         3         4         3         6         6         .         .
    ## G5         1         4        10         5         2         6         8
    ## G6         5         3         6         4         7         5         4

``` r

integrated_mat2[1:6, 1:7]
```

    ## 6 x 7 sparse Matrix of class "dgCMatrix"
    ##    mat1_C1 mat1_C2 mat1_C3 mat1_C4 mat1_C5 mat2_C1 mat2_C2
    ## G1       2       8       4       1       6       .       .
    ## G2       7       3       6       3       5       .       .
    ## G3       5       1       6       4       6       .       .
    ## G4       3       4       3       6       6       .       .
    ## G5       1       4      10       5       2       6       8
    ## G6       5       3       6       4       7       5       4

Key Features:

- Duplicate Resolution: Automatically calls `AggregateDups` to handle
  redundant gene symbols.

- Auto-Naming: Uses argument names (like BatchA) as cell ID prefixes.

**Seurat Integration**

For Seurat objects, `SCIntegrate` automates the standard workflow via
the pipeline parameter.

The Pipeline String: Each letter in the pipeline argument triggers a
specific Seurat command:

| Code  | Function               | Description                         |
|-------|------------------------|-------------------------------------|
| **o** | `CreateSeuratObject`   | (do not use it in this function)    |
| **n** | `NormalizeData`        | Standard normalization.             |
| **s** | `ScaleData`            | Scales data for PCA.                |
| **v** | `FindVariableFeatures` | Selects highly variable genes.      |
| **p** | `RunPCA`               | Principal Component Analysis.       |
| **e** | `FindNeighbors`        | Computes SNN graph.                 |
| **c** | `FindClusters`         | Louvain algorithm clustering.       |
| **t** | `RunTSNE`              | t-SNE reduction.                    |
| **u** | `RunUMAP`              | UMAP reduction.                     |
| **r** | `SCTransform`          | **SCT workflow.** Replaces n, s, v. |

If more are needed, see [A Guide for Extend Single Cell Preprocess
Pipeline](https://wanglabcsu.github.io/SigBridgeR/articles/Extending_Seurat.html)
(deprecated).

Example Usage

``` r

integrated <- SCIntegrate(
  obj1,
  obj2, # -> merge.Seurat
  pipeline = "nsvpi",
  method = Seurat::RPCAIntegration, # Change integration method
  dims = 1:30, # Passed to RunPCA and IntegrateLayers
  k.weight = 50 # Passed to IntegrateLayers
)
```

An example using mock data

``` r

mat1 <- matrix(
  rpois(1000, 5),
  nrow = 20,
  dimnames = list(paste0("G", 1:20), paste0("C", 1:50))
)
mat2 <- matrix(
  rpois(1200, 6),
  nrow = 20,
  dimnames = list(paste0("G", 11:30), paste0("C", 1:60))
)

seu1 <- Seurat::CreateSeuratObject(mat1)
```

    ## Warning: Data is of class matrix. Coercing to dgCMatrix.

``` r

seu2 <- Seurat::CreateSeuratObject(mat2)
```

    ## Warning: Data is of class matrix. Coercing to dgCMatrix.

``` r

integrated_seu <- SCIntegrate(
  Seurat1 = seu1, # Add prefixes
  Seurat2 = seu2,
  method = Seurat::CCAIntegration,
  pipeline = "nsvpi",
  dims = 1:10,
  k.weight = 40
)
```

    ## Start merging Seurat objects

    ## Normalizing layer: counts.1

    ## Normalizing layer: counts.2

    ## Centering and scaling data matrix

    ## Finding variable features for layer counts.1

    ## Finding variable features for layer counts.2

    ## Warning: Requested number is larger than the number of available items (10).
    ## Setting to 10.

    ## Warning: Requested number is larger than the number of available items (10).
    ## Setting to 10.
    ## Warning: Requested number is larger than the number of available items (10).
    ## Setting to 10.
    ## Warning: Requested number is larger than the number of available items (10).
    ## Setting to 10.
    ## Warning: Requested number is larger than the number of available items (10).
    ## Setting to 10.

    ## PC_ 1 
    ## Positive:  G20, G14, G17, G15, G19 
    ## Negative:  G16, G18, G12, G11, G13 
    ## PC_ 2 
    ## Positive:  G17, G11, G19, G15, G13 
    ## Negative:  G18, G14, G20, G12, G16 
    ## PC_ 3 
    ## Positive:  G11, G16, G13, G15, G14 
    ## Negative:  G12, G17, G19, G18, G20 
    ## PC_ 4 
    ## Positive:  G12, G11, G15, G20, G13 
    ## Negative:  G19, G18, G16, G14, G17 
    ## PC_ 5 
    ## Positive:  G15, G16, G17, G18, G11 
    ## Negative:  G13, G12, G19, G20, G14

    ## Finding all pairwise anchors

    ## Warning in RunCCA.Seurat(object1 = object.1, object2 = object.2, assay1 =
    ## "ToIntegrate", : Fewer than 50 features used as input for CCA.

    ## Running CCA

    ## Merging objects

    ## Finding neighborhoods

    ## Finding anchors

    ##  Found 213 anchors

    ## Merging dataset 1 into 2

    ## Extracting anchors for merged samples

    ## Finding integration vectors

    ## Finding integration vector weights

    ## Integrating data

    ## ✔ Integration finished

``` r

integrated_seu
```

    ## An object of class Seurat 
    ## 30 features across 110 samples within 1 assay 
    ## Active assay: RNA (30 features, 10 variable features)
    ##  5 layers present: counts.1, counts.2, data.1, data.2, scale.data
    ##  2 dimensional reductions calculated: pca, integrated.dr

#### 1.2 Bulk expression data

Here are some methods for processing bulk RNA-seq gene expression data
matrices.

##### 1.2.1 Convert bulk RNA-seq raw counts to TPM

``` r

set.seed(123)
ngenes <- 1000
nsamples <- 1000

gene_names <- paste0("GENE", seq_len(ngenes))
sample_names <- paste0("Sample", seq_len(nsamples))

counts <- matrix(
  rnbinom(ngenes * nsamples, mu = 100, size = 10),
  nrow = ngenes,
  ncol = nsamples,
  dimnames = list(gene_names, sample_names)
)
gene_length <- sample(
  500:5000,
  ngenes,
  replace = TRUE
)
names(gene_length) <- gene_names
```

``` r

tpm <- CountsToTPM(counts, gene_length)
tpm[1:5, 1:5]
```

    ##        Sample1   Sample2   Sample3   Sample4   Sample5
    ## GENE1 602.7856 1078.0287  687.1092  514.0647  647.0981
    ## GENE2 224.4745  335.5068  437.5985  556.4223  384.2567
    ## GENE3 760.3323  262.8887  411.2022  511.3512  255.3896
    ## GENE4 537.6406 1100.3214 1685.6017 1223.9918 1116.6870
    ## GENE5 609.4410  452.0651  437.2982  840.6543  329.2632

##### 1.2.2 Gene Symbol Conversion

`SymbolConvert` performs a straightforward task: converting common gene
identifiers (e.g., Ensemble IDs, Entrez) to standardized gene symbols by
using [IDConverter](https://github.com/ShixiangWang/IDConverter)
package.

``` r

# genes * samples
bulk <- read.csv(
  "path_to_your_file.csv",
  header = TRUE,
  row.names = 1
)
```

``` r

bulk <- SymbolConvert(bulk)
```

You can also use other package like `org.Hs.eg.db` or `SeuratExtend` for
gene symbol matching if you prefer not to use `SymbolConvert`’s built-in
`IDConverter`.

``` r

library(org.Hs.eg.db)

bulk <- read.csv(
  "path_to_your_file.csv",
  header = TRUE,
  row.names = 1
)

ensembl_ids <- sub("\\..*", "", rownames(bulk))
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

rownames(bulk) <- gene_symbols
```

``` r

library(SeuratExtend)

bulk <- read.csv(
  "path_to_your_file.csv",
  header = TRUE,
  row.names = 1
)

ensembl_ids <- sub("\\..*", "", rownames(bulk))
gene_symbols <- EnsemblToGenesymbol(
  org.Hs.eg.db,
  spe = "human",
  keep.seq = TRUE
)

rownames(bulk) <- gene_symbols
```

#### 1.3 Phenotype Data

Basically you can just use your phenotype data directly. If you are
confused about the structure
[`Screen()`](https://wanglabcsu.github.io/sigbridger/reference/Screen.md)
requires, please refer to [Section 4 Example](#example).

We provide some functions helping formatting and checking your phenotype
data.

**Checking NA and report the position**:

``` r

# `max_print`: output how many NA location messages at one time in console if NA exists
CheckNA(your_phenotype_data, max_print = 5L)
```

``` r

mat <- matrix(
  c(NA, 1, 1, NA),
  2,
  dimnames = list(c("Gene1", "Gene2"), c("Sample1", "Sample2"))
)
CheckNA(mat)
```

    ## ! Found 2 NA values in data

    ## First 2 positions:

    ## Row 1 ("Gene1"), Col 1 ("Sample1")

    ## Row 2 ("Gene2"), Col 2 ("Sample2")

``` r

mtcars$wt[5] <- NA
CheckNA(mtcars)
```

    ## ! Found 1 NA value in data

    ## First 1 position:

    ## Row 5 ("Hornet Sportabout"), Col 6 ("wt")

**In-place data transformation**

Use it just in tidyverse-like style. The difference between `PhenoMap`
and the dplyr workflow is that `PhenoMap` is slightly faster when
processing two-dimensional `data.frame` data, and it preserves the
`names` attribute when processing one-dimensional vectors. Apart from
these two points, there is no other difference - it is simply a
syntactic sugar we provide.

``` r

# ? if a vector
v <- 1:2000

v2 <- PhenoMap(v, v < 1000 ~ "0", v > 1000 ~ "1", .default = "here_is_1k")
table(v2)
```

    ## v2
    ##          0          1 here_is_1k 
    ##        999       1000          1

``` r

v3 <- PhenoMap(v, v < 100 ~ 0, v < 1000 ~ 1, v > 1000 ~ 2, .default = NA_real_)
table(v3)
```

    ## v3
    ##    0    1    2 
    ##   99  900 1000

``` r

# ? if a data
d <- mtcars
d2 <- PhenoMap(d, mpg > 15 ~ 1, mpg <= 15 ~ 0)
table(d2$mpg)
```

    ## 
    ##  0  1 
    ##  6 26

------------------------------------------------------------------------

### 2. Screen Cells Associated with Phenotype

The function **`Screen`** provide 12 different options for screening
cells associated with phenotype, These 12 algorithms come from the
repositories mentioned in [6. References](#id_6-references), and you can
choose one of them to screen your cells.

Key parameters for `Screen`:

- `matched_bulk`: A data frame of bulk expression data after
  intersecting samples. Make sure the rownames of `matched_bulk` is
  identical to `phenotype`.
- `sc_data`: A Seurat object after preprocessing, you can use the output
  of `Preprocess` function or your own preprocessed Seurat object.
- `phenotype`: A data frame or named vacor of phenotype data after
  intersecting samples. See [5. Example](#id_5-example) for more
  details.
- `label_type`: A character value specifying the filtering labels are
  stored in the `Seurat_object@misc` . Default: `NULL`, meaning the name
  of metho will be used.
- `phenotype_class`: A character value specifying the phenotype data
  type, i.e. `"binary"`, `"survival"` or `"continuous"`.
- `screen_method`: A character value specifying the screening method,
  i.e. “Scissor”, “scPAS”, “scAB”, “scPP”, “DEGAS”, “DEGASv2”, “LP_SGL”,
  “PIPET”, “SIDISH”, “SCIPAC”, “TiRank”, or “PhenoMapR”. Use
  `AVAILABLE_SCREEN_METHODS` to list all registered methods.
- `...`: Other parameters for the screening methods, see below

#### 2.1 (Option A) Scissor Screening

Parameters pass to `...` when using `Scissor` method:

- `save_path`: A character value specifying the path to save
  intermediate data, you can also set `NULL` to suppress the saving of
  intermediate files. Default:
  `Scissor_res/<cache_stamp>/Scissor_inputs.RData`
- `load_path`: A character value specifying the path to load
  intermediate data, e.g. `"Scissor_res/<cache_stamp>`. Default: `NULL`
- `alpha`: Parameter used to balance the effect of the l1 norm and the
  network-based penalties. It can be a number or a searching vector. If
  alpha = NULL, a default searching vector is used. The range of alpha
  is in `[0,1]`. A larger alpha lays more emphasis on the l1 norm.
- `cutoff`: Cutoff for the percentage of the Scissor selected cells in
  total cells. This parameter is used to restrict the number of the
  Scissor selected cells. A cutoff less than 50% (default 20%) is
  recommended depending on the input data. Only used when
  `alpha = NULL`.
- `reliability_test`: A logical value specifying whether to perform
  reliability test. Default: `FALSE`
- `reliability_test.n`: Permutation times (default: 10)
- `reliability_test.nfold`: The fold number in cross-validation
  (default: 10)
- `cell_evaluation`: A logical value specifying whether to perform cell
  evaluation. Default: `FALSE`
- `cell_evaluation.benchmark_data`: Path to benchmark data (RData file).
- `cell_evaluation.FDR`: FDR threshold for cell evaluation (default:
  0.05).
- `cell_evaluation.bootstrap_n`: Number of bootstrap iterations for cell
  evaluation (default: 100).

**Usage**:

``` r

scissor_result <- Screen(
  matched_bulk = matched_bulk,
  sc_data = sc_dataset, # A Seurat object after preprocessing
  phenotype = matched_phenotype_data,
  label_type = "TP53", # The filtering labels are stored in the `@misc`, you can change it to your own label
  phenotype_class = "binary",
  screen_method = c("Scissor"),
  save_cache = "Tmp" # Intermediate data
)
```

You can use the intermediate data for repeated runs. This is an inherent
feature of the `Scissor`.

``` r

scissor_result <- Screen(
  sc_data = sc_dataset,
  label_type = "TP53",
  phenotype_class = "binary",
  screen_method = c("Scissor"),
  load_cache = "Tmp" # Intermediate data
)
```

If only the parameters `alpha` and `cutoff` are adjusted, this method
can also be applied.

``` r

# When `alpha = NULL`, an alpha iteration will continue until phenotype-associated cells are screened out or no cells are screened out even after exceeding the `cutoff`.
scissor_result <- Screen(
  sc_data = sc_dataset,
  label_type = "TP53",
  phenotype_class = "binary",
  screen_method = c("Scissor"),
  load_cache = "Tmp", # Intermediate data
  alpha = NULL,
  cutoff = 0.2
)
```

**Returning structure**: A list containing:

- `scRNA_data`: A Seurat object after screening
- `scissor_result`: The result of Scissor screening, including the
  parameters used
- `reliability_test`: Reliability test results
- `cell_evaluation`: Cell evaluation results

**Cell level Evaluation & Reliability Test**:

You can use `cell_evalutaion = TRUE` and `reliability_test = TRUE` to
obtain some supporting information for each Scissor selected cell.
First, prepare a benchmark dataset yourself for cell evalutaion.

``` r

scissor_result <- Screen(
  sc_data = sc_dataset,
  label_type = "TP53",
  phenotype_class = "binary",
  screen_method = c("Scissor"),
  load_cache = "Tmp", # Intermediate data
  reliability_test = TRUE,
  cell_evaluation = TRUE,
  cell_evaluation.benchmark_data = "path_to_benchmark_data.RData",
  alpha = NULL,
  cutoff = 0.05
)
```

helpful documentation:

[Scissor-Cell Level
Evaluations](https://sunduanchen.github.io/Scissor/vignettes/Scissor_Tutorial.html#cell-level-evaluations)

#### 2.2 (Option B) scPAS Screening

Parameters pass to `...` when using `scPAS` method (basically adapted
from the `scPAS`’s documentation):

- Parameters passed to `scPAS::scPAS()`

  These parameters directly interface with the core `scPAS`() function
  from the original package:

  - `assay`: Name of Assay to get.
  - `imputation`: Logical. imputation or not.
  - `imputation_method`: Character. Name of alternative method for
    imputation.
  - `nfeature`: Numeric. The Number of features to select as top
    variable features in `sc_data`. Top variable features will be used
    to intersect with the features of `matched_bulk`. Default is NULL
    and all features will be used.
  - `alpha`: Numeric. Parameter used to balance the effect of the l1
    norm and the network-based penalties. It can be a number or a
    searching vector. If `alpha = NULL`, a default searching vector is
    used. The range of alpha is in `[0,1]`. A larger alpha lays more
    emphasis on the l1 norm.
  - `cutoff`: Numeric. Cutoff for the percentage of the scPAS selected
    cells in total cells when `alpha = NULL`. This parameter is used to
    restrict the number of the scPAS selected cells. A cutoff less than
    50% (default 20%) is recommended depending on the input data.
  - `network_class`: The source of feature-feature similarity network.
    By default this is set to sc and the other one is bulk.
  - `FDR_threshold`: Numeric. FDR value threshold for identifying
    phenotype-associated cells (default: 0.05)
  - `independent`: Logical. The background distribution of risk scores
    is constructed independently of each cell. (default: TRUE)
  - `permutation_times`: Number of permutations to perform (default:
    2000)

**usage**:

``` r

scpas_result <- Screen(
  matched_bulk = matched_bulk,
  sc_data = A_Seurat_object,
  phenotype = phenotype,
  label_type = "TP53", # The filtering labels are stored in the `@misc`
  screen_method = "scPAS",
  phenotype_class = "binary",
  assay = 'RNA',
  imputation = FALSE,
  imputation_method = c("KNN", "ALRA"),
  nfeature = 3000L,
  alpha = c(0.01, NULL),
  cutoff = 0.2,
  network_class = c("SC", "bulk"),
  permutation_times = 2000L,
  FDR_threshold = 0.05,
  independent = TRUE
)
```

**returning structure**: A list containing:

- `scRNA_data`: A Seurat object after screening
- `stats`: A data.frame the significance of scPAS screening results
- `para`: A list containing the parameters used in scPAS screening

#### 2.3 (Option C) scAB Screening

Parameters pass to `...` when using `scAB` method (basically adapted
from the `scAB`’s documentation):

- `alpha`: Coefficient of phenotype regularization, default is `0.005`.
  When specified `NULL`, a default searching vector is used. A custom
  numeric vector is also supported.
- `alpha_2`: Coefficient of cell-cell similarity regularization, default
  is `0.005`. When specified `NULL`, a default searching vector is used.
  A custom numeric vector is also supported
- `maxiter`: Maximum number of iterations, default is `2000`
- `tred`: Threshold for early stopping, default is `2`

**usage**:

``` r

scab_result <- Screen(
  matched_bulk = your_matched_bulk,
  sc_data = A_Seurat_object,
  phenotype = your_matched_phenotype,
  label_type = "TP53", # The filtering labels are stored in the `@misc`
  screen_method = "scAB",
  phenotype_class = "binary",
  alpha = c(0.005, NULL),
  alpha_2 = c(0.005, NULL),
  maxiter = 2000L,
  tred = 2L
)
```

**note**:

When both `alpha` and `alpha_2` are specified as `NULL` or as search
vectors, the total number of searches equals the product of their
lengths, which may lead to long runtimes. In such cases, we recommend
enabling parallel computation.

**returning structure**: A list containing:

- `scRNA_data`: A Seurat object after screening
- `scAB_result`: A list with the submatrix and loss value

#### 2.4 (Option D) scPP Screening

Parameters pass to `...` when using `scPP` method :

- `ref_group`: The reference group for the binary analysis, default is
  `0`
- `Log2FC_cutoff`: The cutoff for the log2 fold change of the binary
  analysis, default is `0.585`
- `estimate_cutoff`: Effect size threshold for continuous traits,
  default is `0.2`
- `probs`: Quantile cutoff in (0, 0.5) for cell classification, default
  is `0.2`. When specified `NULL`, a default searching vector is used. A
  custom searching vector is also supported.

**usage**:

``` r

# This will take several hours
scpp_result <- Screen(
  matched_bulk = your_matched_bulk,
  sc_data = A_Seurat_object,
  phenotype = your_matched_phenotype,
  label_type = "TP53", # The filtering labels are stored in the `@misc`
  screen_method = "scpp",
  phenotype_class = "binary",
  ref_group = 0,
  Log2FC_cutoff = 0.585,
  estimate_cutoff = 0.2,
  probs = c(0.2, NULL)
)
```

**returning structure**: A list containing:

- `scRNA_data`: A Seurat object after screening
- `gene_list`: A list containing positive genes and negative genes

#### 2.5 (Option E) DEGAS Screening

Parameters pass to `...` when using `DEGAS` method

- `sc_data.pheno_colname`: The column name of the phenotype in the
  `sc_data@meta.data` slot, used to specify the phenotype for more
  accurate screening. Default is `NULL`.
- `tmp_dir`: The directory for storing the intermediate files. Default
  is `NULL`, a directory named `tmp` will be created.
- `env_params` : A list of parameters for the environment, default is
  [`list()`](https://rdrr.io/r/base/list.html). Use
  [`?DoDEGAS`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md)
  to see the details.
- `degas_params`: A list of parameters for the DEGAS algorithm, default
  is [`list()`](https://rdrr.io/r/base/list.html). Use
  [`?DoDEGAS`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md)
  to see the details.
- `normality_test_method`: Method for normality testing:
  `"jarque-bera", "d'agostino", or "kolmogorov-smirnov"`, default is
  “jarque-bera”.

``` r

reticulate::use_condaenv("r-reticulate-degas")

degas_result <- Screen(
  matched_bulk = your_matched_bulk,
  sc_data = A_Seurat_object,
  phenotype = your_matched_phenotype,
  label_type = "TP53", # The labels are stored in the `@misc` and are used to identify the screening results.
  screen_method = "DEGAS",
  # The type of phenotype
  phenotype_class = c("binary", "continuous", "survival"),
  # DEGAS parameters
  degas_params = list(
    # DEGAS.model_type will be automatically determined by the `phenotype_class`
    DEGAS.model_type = c(
      "BlankClass", # only bulk level phenotype specified
      "ClassBlank", # only single cell level phenotype specified
      "ClassClass", # when both single cell level phenotype and bulk level phenotype specified
      "ClassCox", # when both single cell level phenotype and bulk level survival data specified
      "BlankCox" # only bulk level survival data specified
    ),
    # DEGAS.architecture will be `DenseNet` in default
    DEGAS.architecture = c(
      "DenseNet", # a dense net network
      "Standard" # a feed forward network
    ),
    # The path to save intermediate data
    path.data = '',
    path.result = '',
    # The python executable path of the conda environment, auto detected in default
    DEGAS.pyloc = NULL,
    # Some functions will be called in the DEGAS algorithm
    DEGAS.toolsPath = paste0(.libPaths()[1], "/DEGAS/tools/"),
    # Screening parameters
    DEGAS.ff_depth = 3,
    DEGAS.bag_depth = 5,
    DEGAS.train_steps = 2000,
    DEGAS.scbatch_sz = 200,
    DEGAS.patbatch_sz = 50,
    DEGAS.hidden_feats = 50,
    DEGAS.do_prc = 0.5,
    DEGAS.lambda1 = 3.0,
    DEGAS.lambda2 = 3.0,
    DEGAS.lambda3 = 3.0,
    DEGAS.seed = 2
  ),
  # default: "jarque-bera".
  normality_test_method = c(
    "jarque-bera",
    "d'agostino",
    "kolmogorov-smirnov"
  )
)
```

After update to v4.0.0. `Screen` does not automatically creat a python
env for DEGAS, you need to create a python env for DEGAS before using
`Screen`. So as to comply the R package’s principle.

To obtain a template python env configuration file for DEGAS, you can
use:

``` r

env_config <- system.file(
  "conda/DEGAS_environment.yml",
  package = "SigBridgeR"
) |>
  readLines()
```

To obtain the default parameters for DEGAS, you can use:

``` r

degas_param <- SigBridgeR:::DEGASParamSet(list())

env_param <- SigBridgeR:::DEGASEnvSet(list())
```

You can use `ListPyEnvs()` to list all the python environments in your
system, including virtual environments. Both Windows and Unix-like
systems are supported. More information can be found in
[Other_Function_Details](https://wanglabcsu.github.io/SigBridgeR/articles/Other_Function_Details.html)

``` r

# * On Unix-like system it goes like this
ListPyEnv()
#                 name                                                   python  type
# 1               base                         /home/user/miniconda3/bin/python conda
# 2 r-reticulate-degas /home/user/miniconda3/envs/r-reticulate-degas/bin/python conda
# 3               test                  /home/user/.virtualenvs/test/bin/python  venv
```

Please note that the environmental dependencies required for DEGAS to
run are quite stringent, and conflicts are highly likely to occur.

**returning structure**: A list containing:

- `scRNA_data`: A Seurat object after screening
- `model`: A model trained using single-cell RNA expression matrix,
  tissue bulk RNA sequencing expression matrix, and phenotypic data.
- `DEGAS_prediction`: Using the model to conduct prediction for each
  cell, resulting in a data.frame where each phenotype has a predicted
  probability score.

#### 2.6 (Option F) LP_SGL Screening

Parameters pass to `...` when using `LP_SGL` method

- `resolution`: Resolution parameter for Leiden clustering (default:
  `0.6`)
- `alpha`: Alpha parameter for SGL balancing L1 and L2 penalties
  (default: `0.5`)
- `nfold`: Number of folds for cross-validation (default: `5`)
- `dge_analysis`: List controlling differential expression analysis:
  - `run`: Whether to run DEG analysis (default: `FALSE`)
  - `logFC_threshold`: Log fold change threshold (default: `1`)
  - `pval_threshold`: P-value threshold (default: `0.05`)

``` r

lpsgl_result <- Screen(
  matched_bulk = your_matched_bulk,
  sc_data = A_Seurat_object,
  phenotype = your_matched_phenotype,
  label_type = "TP53",
  resolution = 0.6,
  alpha = 0.5,
  nfold = 5,
  dge_analysis = list(
    run = FALSE, # whether to run DEG analysis
    logFC_threshold = 1,
    pval_threshold = 0.05
  ),
  ... # Additional parameters like verbose, seed
)
```

**returning structure**: A list containing:

- `scRNA_data`: A Seurat object after screening
- `sgl_fit` : Fitted SGL model object
- `cvfit` : Cross-validation results
- `dge_res` : Differential expression results if requested (NULL
  otherwise)

#### 2.7 (Option G) PIPET Screening

Parameters pass to `...` when using `PIPET` method

- **Phenotype adaptation parameters**:
  - `discretize_method`: Character: “kmeans”/“quantile”/“custom”
    (default: “kmeans”)
  - `cutoff`: Numeric vector for custom discretization when
    `discretize_method` is “custom” (default: NULL)
- **Marker generation parameters**:
  - `log2FC`: Numeric: log2FC threshold (default: 1)
  - `p.adjust`: Numeric: adjusted p-value threshold (default: 0.05)
- **Single-cell annotation parameters**:
  - `distance`: Character: “cosine”/“pearson”/“spearman” (default:
    “cosine”)
  - `nPerm`: Integer: permutation times for statistical test (default:
    1000L)

**Usage**:

``` r

pipet_result <- Screen(
    matched_bulk = matched_bulk,
    sc_data = sc_dataset, # A Seurat object after preprocessing
    phenotype = matched_phenotype_data,
    label_type = "TP53", # The filtering labels are stored in the `Seurat_object@misc`
    phenotype_class = "binary", # `survival` is not supported
    screen_method = "PIPET",
    # PIPET specific parameters
    lg2FC = 1,
    p.adjust = 0.05,
    distance = "cosine",
    nPerm = 1000,
    parallel = FALSE # Whether to use parallel computing, before using this parameter, please make sure that future::plan() has been called
)
```

**Returning structure**: A list containing:

- `scRNA_data`: A Seurat object after screening (with PIPET annotations
  in meta.data)
- `markers`: Phenotype-specific marker genes

#### 2.8 (Option H) SIDISH Screening

Parameters passed to `...` when using `SIDISH` method:

- `verbose`: Logical. Whether to print verbose output during execution.
  Default is `TRUE`.
- `label_type`: Character specifying the phenotype label type stored in
  `sc_data@misc`. Default is `"SIDISH"`.
- `assay`: Seurat assay name to use for expression data. Default is
  `"RNA"`.
- `sidish_param`: List of parameters for SIDISH algorithm. See below for
  details.
- `env_params`: List of parameters for environment setup. See below for
  details.

``` r

sidish_result <- Screen(
  matched_bulk = your_matched_bulk,
  sc_data = A_Seurat_object,
  phenotype = your_survival_phenotype, # Must contain "time" and "status" columns
  label_type = "SIDISH", # Labels stored in `@misc` slot to identify screening results
  screen_method = "SIDISH",
  phenotype_class = "survival", # Currently only survival phenotype is supported

  # SIDISH algorithm parameters
  sidish_param = list(
    # Preprocessing parameters
    patient_id = "Sample",
    celltype_name = "celltype_major",
    processed = TRUE,
    n_genes_by_counts = 5000L,
    pct_counts_mt = 10L,
    batch_correction = FALSE,
    survival_ = "time", # Column name for survival time
    status = "status", # Column name for event status

    # Execution environment
    device = "cuda", # "cuda" for GPU acceleration, "cpu" for CPU-only
    use_spatial_graph = FALSE,
    k_neighbors = NULL,

    # Phase 1: VAE training parameters
    phase1_epochs = 225L,
    phase1_i_epochs = 20L,
    phase1_latent_size = 32L,
    phase1_layer_dims = c(512L, 128L),
    phase1_batch_size = 256L,
    phase1_optimizer = "Adam",
    phase1_lr = 1e-4,
    phase1_lr_3 = 1e-4,
    phase1_dropout = 0L,
    phase1_type = "Dense", # "Dense" or "Normal"

    # Phase 2: Deep Cox training parameters
    phase2_epochs = 500L,
    phase2_hidden = 128L,
    phase2_lr = 1e-4,
    phase2_dropout = 0L,
    phase2_test_size = 0.2,
    phase2_batch_size_bulk = 256L,

    # Training & risk definition
    train_iterations = 5L,
    train_percentile = 0.95,
    train_steepness = 30L,
    train_path = "./SIDISH_res/", # Path to save intermediate data
    train_num_workers = 0L,
    train_distribution_fit = "fitted" # "fitted" or "default"
  ),

  verbose = TRUE
)
```

**Note**: SIDISH currently **only supports survival phenotypes**
(`phenotype_class = "survival"`). The phenotype data frame must contain
columns named `"time"` (survival time) and `"status"` (event indicator).

To obtain the default parameters for SIDISH, use
`SigBridgeR::SIDISHParamSet()`

``` r

sidish_default_params <- SigBridgeR:::SIDISHParamSet(list())
```

You can use `ListPyEnvs()` to list all Python environments available on
your system (both conda and virtualenv):

``` r

ListPyEnv()
#                 name                                                   python  type
# 1               base                         /home/user/miniconda3/bin/python conda
# 2 r-reticulate-sidish-nvidia /home/user/miniconda3/envs/r-reticulate-sidish-nvidia/bin/python conda
# 3 r-reticulate-sidish-cpu    /home/user/miniconda3/envs/r-reticulate-sidish-cpu/bin/python conda
```

**Return structure**: A named list containing:

- `scRNA_data`: A Seurat object with SIDISH screening results integrated
  into the `@misc` and `@meta.data` slots under the specified
  `label_type`. The object contains cell-level risk scores and
  survival-related annotations generated by the SIDISH algorithm.

#### 2.9 (Option I) SCIPAC Screening

Parameters pass to `...` when using `SCIPAC` method

- `hvg`: Integer. Number of highly variable genes to use for
  preprocessing. Default is `1000L`.
- `do_pca_sc`: Logical. Whether to perform PCA on single-cell data and
  apply the rotation matrix to bulk data; if FALSE, PCA is performed on
  bulk data and applied to single-cell data. Default is `FALSE`.
- `n_pc`: Integer. Number of principal components to use. Default is
  `60L`.
- `sc_batch_col`: Character or vector. Batch variable for single-cell
  data. Default is `NULL`.
- `resolution`: Integer. Clustering resolution for cell type
  identification. Default is `2L`.
- `ela_net_alpha`: Numeric. Elastic net mixing parameter (0 = ridge, 1 =
  lasso). Default is `0.4`.
- `bt_size`: Integer. Bootstrap sample size for stability assessment.
  Default is `50L`.
- `ncore`: Integer. Number of CPU cores for parallel computation.
  Default is `7L`.
- `ci_alpha`: Numeric. Significance level for confidence intervals.
  Default is `0.05`.
- `nfold`: Integer. Number of folds for cross-validation for regression
  models. Default is `10L`.
- `...`: Additional arguments. Supports `assay` (character), `verbose`
  (logical), and `seed` (integer).

``` r

scipac_result <- Screen(
  matched_bulk = your_matched_bulk,
  sc_data = A_Seurat_object,
  phenotype = your_matched_phenotype,
  label_type = "TP53",
  screen_method = "SCIPAC",
  phenotype_class = "binary",
  hvg = 1000L,
  do_pca_sc = FALSE,
  n_pc = 60L,
  sc_batch_col = NULL,
  resolution = 2L,
  ela_net_alpha = 0.4,
  bt_size = 50L,
  ncore = 7L,
  ci_alpha = 0.05,
  nfold = 10L,
)
```

**returning structure**: A list containing:

- `scRNA_data`: A Seurat object after screening
- `pca_res` : PCA rotation results
- `cluster_res` : Clustering results of single-cell data

#### 2.10 (Option J) DEGASv2 Screening

Parameters pass to `...` when using `DEGASv2` method (basically adapted
from the `DEGASv2`’s documentation):

**Feature selection (shared by bulk and single-cell data)**

- `bulk_hvg`: Logical. Whether to select highly variable genes from the
  bulk data. Default is `TRUE`.
- `bulk_de`: Logical. Whether to perform differential expression
  analysis on the bulk data. Default is `TRUE`.
- `sc_de`: Logical. Whether to perform differential expression analysis
  on the single-cell data. Default is `TRUE`.
- `add_genes`: Character vector or `NULL`. Additional genes forced into
  the feature set. Default is `NULL`.
- `n_hvg`: Integer. Number of highly variable genes to select from bulk
  data. Default is `250L`.
- `n_bulk_de`: Integer. Number of bulk DE genes to select. Default is
  `250L`.
- `n_sc_de`: Integer. Number of single-cell DE genes to select. Default
  is `200L`.
- `padj.thresh`: Numeric. Adjusted p-value threshold for DE gene
  selection. Default is `0.05`.
- `only.pos`: Logical. Whether to keep only positively enriched markers.
  Default is `FALSE`.
- `min.pct`: Numeric. Minimum fraction of cells expressing a gene in the
  DE analysis. Default is `0.25`.
- `logfc.threshold`: Numeric. Log fold-change threshold for the DE
  analysis. Default is `0.25`.

**Cell labelling**

- `top_fraction_pos`: Numeric. Proportion of cells labelled as
  `"Positive"` according to the predicted hazard (the top quantile).
  Default is `0.2`.
- `sclab`: Character or `NULL`. Metadata column of `sc_data` providing
  supervised cell labels. When supplied, both single-cell and bulk
  labels are used (`ClassClass` / `ClassCox` mode). Default is `NULL`.
- `n_st_classes`: Integer. Number of single-cell classes. Default is
  `length(unique(sc_data$seurat_clusters))`.

**Model training**

- `loss_type`: Character. Loss function, one of `"cross_entropy"`,
  `"log_neg"`, `"rank_loss"`. Default is `"cross_entropy"`.
- `transfer_type`: Character. Domain adaptation method, one of
  `"Wasserstein"`, `"MMD"`. Default is `"Wasserstein"`.
- `model_save_dir`: Character. Directory used to store the trained
  model. Default is `"DEGASv2_res"`.
- `lambda1`, `lambda2`, `lambda3`: Numeric. Weights of the
  classification, domain-adaptation and reconstruction losses. Defaults
  are `1L`, `3L` and `3L`.
- `tot_seeds`: Integer. Number of random seeds to train. Default is
  `10L`.
- `tot_iters`: Integer. Number of training iterations. Default is
  `300L`.
- `extract_embs`: Logical. Whether to extract cell embeddings. Default
  is `FALSE`.
- `random_feat`: Logical. Whether to use random features. Default is
  `FALSE`.
- `random_perc`: Numeric. Proportion of random features when
  `random_feat = TRUE`. Default is `0.8`.
- `early_stopping`: Logical. Whether to enable early stopping. Default
  is `FALSE`.

**Usage**:

``` r

degasv2_result <- Screen(
  matched_bulk = your_matched_bulk, # raw counts matrix is preferred
  sc_data = A_Seurat_object,
  phenotype = your_matched_phenotype,
  label_type = "TP53",
  screen_method = "DEGASv2",
  phenotype_class = "binary", # "binary", "survival" or "continuous"
  # Cell labelling
  top_fraction_pos = 0.2,
  # Optional supervised single-cell labels
  sclab = NULL,
  # Feature selection
  bulk_hvg = TRUE,
  bulk_de = TRUE,
  sc_de = TRUE,
  add_genes = NULL,
  n_hvg = 250L,
  n_bulk_de = 250L,
  n_sc_de = 200L,
  padj.thresh = 0.05,
  # Training
  loss_type = "cross_entropy",
  transfer_type = "Wasserstein",
  model_save_dir = "DEGASv2_res",
  lambda1 = 1L,
  lambda2 = 3L,
  lambda3 = 3L,
  tot_seeds = 10L,
  tot_iters = 300L,
  early_stopping = FALSE
)
```

``` r

table(degasv2_result$scRNA_data$DEGASv2)
#    Other Positive
```

**note**:

- DEGASv2 is the re-implemented, deep-transfer-learning based successor
  of DEGAS. Its training backend is PyTorch, accessed through
  `reticulate`, but unlike DEGAS you do **not** need to create and
  activate a Python environment manually: the `DEGASv2` package declares
  its Python dependencies with
  [`reticulate::py_require()`](https://rstudio.github.io/reticulate/reference/py_require.html),
  so they are provisioned automatically on first use. The `DEGASv2`
  package itself must be installed (`pak::pak("Exceret/DEGASv2")`,
  handled automatically by `Screen`).
- Bulk raw counts are preferred, since the bulk branch can perform
  differential expression analysis (`bulk_de = TRUE`); if the bulk
  matrix is not a counts matrix, the bulk DE step is skipped with a
  warning.

**returning structure**: A list containing:

- `scRNA_data`: A Seurat object after screening. The screening result is
  stored in the `DEGASv2` column (`"Positive"` / `"Other"`), together
  with the predicted hazard (`DEGASv2_hazard`) and the threshold used
  for labelling (`DEGASv2_pos_thresh`). The parameters used are recorded
  in `scRNA_data@misc$DEGASv2`.

#### 2.11 (Option K) TiRank Screening

TiRank integrates bulk and single-cell data through a rank-based deep
learning framework. Parameters pass to `...` when using `TiRank` method:

- `tirank_params`: List of TiRank algorithm parameters. See below for
  details. Default is the full default list.
- `save_path`: (Soft-deprecated) Directory for intermediate and final
  results. The final prediction table is read from
  `<save_path>/3_Analysis/spot_predict_score.csv`, so keep it aligned
  with the cache directory when caching is used. Default is
  `"./TiRank_res"`.
- `assay`: Character. Assay of `sc_data` to use. Default is `"RNA"`.
- `load_cache` / `save_cache`: Character. Cache directory used to
  load/save intermediate results (preferred over `save_path`), see
  [`CacheSetHere()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSetHere.md).
  Default is `NULL` / `save_path`.
- `verbose`: Logical. Whether to print progress messages. Default is
  inherited from `getFuncOption("verbose")`.
- `seed`: Integer. Random seed. Default is inherited from
  `getFuncOption("seed")`.

Parameters in `tirank_params`:

**Data preprocessing**

- `validation_proportion`: Numeric. Proportion of bulk data held out for
  validation. Default is `0.15`.
- `sampling_thresh`: Numeric. Threshold for the resampling strategy.
  Default is `0.5`.
- `sampling_mode`: Character. Resampling method, one of `"smote"`,
  `"downsample"`, `"upsample"`, `"tomeklinks"`.
- `top_var_genes`: Integer. Number of top variable genes to select.
  Default is `2000L`.
- `top_gene_pairs`: Integer. Number of top gene pairs used for ranking.
  Default is `1000L`.
- `p_value_threshold`: Numeric. P-value threshold for feature selection.
  Default is `0.05`.
- `max_cutoff`, `min_cutoff`: Numeric. Upper/lower cutoffs for
  correlation filtering. Defaults are `0.8` and `-0.8`.

**Neural network architecture**

- `nhead`: Integer. Number of attention heads. Default is `2L`.
- `nhid1`, `nhid2`: Integer. Sizes of the first and second hidden
  layers. Defaults are `96L` and `8L`.
- `n_output`: Integer. Output dimension of the encoder. Default is
  `32L`.
- `nlayers`: Integer. Number of encoder layers. Default is `3L`.
- `n_pred`: Integer. Number of prediction heads. Default is `2L`.
- `dropout`: Numeric. Dropout rate for regularization. Default is `0.5`.
- `encoder_type`: Character. Encoder architecture, one of `"MLP"`,
  `"Transformer"`, `"DenseNet"`.
- `infer_mode`: Character. Inference mode, `"SC"` (single-cell) or
  `"ST"` (spatial transcriptomics).

**Training**

- `n_trials`: Integer. Number of repeated training trials. Default is
  `5L`.
- `do_reject`: Logical. Whether to apply rejection criteria to uncertain
  predictions. Default is `TRUE`.
- `tolerance`: Numeric. Tolerance threshold for rejection. Default is
  `0.05`.

**Usage**:

``` r

tirank_result <- Screen(
  matched_bulk = your_matched_bulk,
  sc_data = A_Seurat_object,
  phenotype = your_matched_phenotype,
  label_type = "TP53",
  screen_method = "TiRank",
  phenotype_class = "survival", # "binary" (Classification), "survival" (Cox), "continuous" (Regression)
  assay = "RNA",
  tirank_params = list(
    # Data preprocessing
    validation_proportion = 0.15,
    sampling_thresh = 0.5,
    sampling_mode = "smote", # "smote", "downsample", "upsample", "tomeklinks"
    top_var_genes = 2000L,
    top_gene_pairs = 1000L,
    p_value_threshold = 0.05,
    max_cutoff = 0.8,
    min_cutoff = -0.8,
    # Neural network architecture
    nhead = 2L,
    nhid1 = 96L,
    nhid2 = 8L,
    n_output = 32L,
    nlayers = 3L,
    n_pred = 2L,
    dropout = 0.5,
    encoder_type = "MLP", # "MLP", "Transformer", "DenseNet"
    infer_mode = "SC", # "SC" or "ST"
    # Training
    n_trials = 5L,
    do_reject = TRUE,
    tolerance = 0.05
  )
)
```

TiRank supports caching, so the preprocessed data and the trained model
can be reused. Note that the final prediction table is read from
`<save_path>/3_Analysis/spot_predict_score.csv`, so keep `save_path`
aligned with the cache directory being used:

``` r

# First run: train the model and save the cache under `TiRank_res/`
tirank_result <- Screen(
  matched_bulk = your_matched_bulk,
  sc_data = A_Seurat_object,
  phenotype = your_matched_phenotype,
  label_type = "TP53",
  screen_method = "TiRank",
  phenotype_class = "survival",
  save_path = "TiRank_res",
  save_cache = "TiRank_res"
)

# Later runs: reuse the cache directory created in the previous step,
# e.g. `TiRank_res/survival_20250101000000`
tirank_result <- Screen(
  matched_bulk = your_matched_bulk, # not required when the cache is complete
  sc_data = A_Seurat_object,
  phenotype = your_matched_phenotype, # not required when the cache is complete
  label_type = "TP53",
  screen_method = "TiRank",
  phenotype_class = "survival",
  save_path = "TiRank_res/survival_20250101000000",
  load_cache = "TiRank_res/survival_20250101000000"
)
```

``` r

table(tirank_result$scRNA_data$TiRank)
# Negative  Neutral Positive
```

**note**:

- TiRank currently runs on **GPU** (`device = "cuda"` is used
  internally), so a CUDA-capable device is required. The `rTiRank`
  package must be installed (`pak::pak("Exceret/rTIRank")`, handled
  automatically by `Screen`), and its Python dependencies are
  provisioned automatically through
  [`reticulate::py_require()`](https://rstudio.github.io/reticulate/reference/py_require.html).

**returning structure**: A list containing:

- `scRNA_data`: A Seurat object after screening, with `TiRank_Reject`,
  `TiRank_Rank_Score` and `TiRank` (`"Positive"`, `"Neutral"`,
  `"Negative"`) added to `meta.data`, and the parameters stored in
  `scRNA_data@misc$TiRank`.
- `cell_cell_distance`: The computed cell-cell similarity/distance
  matrix.

#### 2.12 (Option L) PhenoMapR Screening

PhenoMapR derives a gene-level reference from the matched bulk
expression and phenotype, then scores every single cell and labels the
two tails of the score distribution. Parameters pass to `...` when using
`PhenoMapR` method:

**Phenotype handling**

- `phenotype_column`: Character or `NULL`. Outcome column name for
  binary/continuous phenotypes. Default is `NULL` (uses
  `"phenotype_value"`).
- `survival_time`, `survival_event`: Character. Column names of the
  survival phenotype data frame. Defaults are `"time"` and `"status"`.
- `binary_positive_reference`: Character. Which level of a binary
  phenotype is treated as the reference-positive direction, `"second"`
  or `"first"`. Default is `"second"`.
- `reference_sign`: Integer. Sign applied to the reference scores, use
  `-1` to mirror the score axis. Default is `1L`.

**Bulk reference construction**

- `normalize`: Logical. Whether to normalize the bulk expression.
  Default is `TRUE`.
- `platform`: Character. Bulk platform, one of `"auto"`, `"rnaseq"`,
  `"microarray"`. Default is `"auto"`.
- `probe_annotation`: Probe-to-gene annotation table for microarray
  data. Default is `NULL`.
- `hugo_species`: Character. Species passed to HGNChelper, `"human"` or
  `"mouse"`. Default is `c("human", "mouse")`.

**Single-cell scoring**

- `assay`, `layer`: Character. Seurat assay and layer used for scoring.
  Defaults are `"RNA"` and `"data"`.
- `z_score_cutoff`: Numeric. Z-score cutoff for limiting paired genes.
  Default is `2L`.
- `score_mode`: Character. Scoring mode, one of `"weighted_sum"`,
  `"activity_adjusted"`.
- `vars_to_regress`: Character vector. Variables regressed out before
  scoring. Default is `c("S.Score", "G2M.Score", "nCount_RNA")`.
- `cancer_type`: Character or `NULL`. Optional cancer-type context.
  Default is `NULL`.
- `pseudobulk`: Logical. Whether to score pseudobulk profiles instead of
  single cells. Default is `FALSE`.
- `group_by`: Character or `NULL`. Grouping variable used for pseudobulk
  aggregation. Default is `NULL`.
- `sample_metadata`: Data frame or `NULL`. Sample-level metadata used
  with `pseudobulk = TRUE`. Default is `NULL`.

**Group definition**

- `percentile`: Numeric. Tail percentile used to define the phenotype
  groups. Default is `0.05`.
- `permutation_n`: Integer or `NULL`. Number of permutations. Default is
  `NULL`.
- `permutation_seed`: Integer or `NULL`. Permutation seed, defaults to
  `seed`.

**Usage**:

``` r

phenomapr_result <- Screen(
  matched_bulk = your_matched_bulk,
  sc_data = A_Seurat_object,
  phenotype = your_matched_phenotype,
  label_type = "PhenoMapR",
  screen_method = "PhenoMapR",
  phenotype_class = "survival", # "binary", "survival" or "continuous"
  # Phenotype handling
  survival_time = "time",
  survival_event = "status",
  # Bulk reference construction
  normalize = TRUE,
  platform = "auto",
  hugo_species = "human",
  # Single-cell scoring
  assay = "RNA",
  layer = "data",
  z_score_cutoff = 2L,
  score_mode = "weighted_sum",
  vars_to_regress = c("S.Score", "G2M.Score", "nCount_RNA"),
  # Group definition
  percentile = 0.05
)
```

``` r

table(phenomapr_result$scRNA_data$PhenoMapR, useNA = "ifany")
# Negative  Neutral Positive
```

**note**:

- PhenoMapR labels the upper tail of the score distribution as
  `"Most Adverse"` and the lower tail as `"Most Favorable"`, which are
  mapped to the package-wide vocabulary as `"Positive"` and
  `"Negative"`; cells in the middle are labelled `"Neutral"`.
- Higher scores follow the phenotype-valued-as-`1` direction when the
  reference is derived with the default orientation. For binary
  phenotypes this means `binary_positive_reference = "second"` regresses
  on the `"1"` level, and for survival traits the Cox z-score is
  positive for hazard-increasing genes. Therefore, when phenotype `1`
  denotes the unfavourable direction, the `"Most Adverse"` tail is the
  one labelled `"Positive"`. If your phenotype is encoded the other way
  round, adjust `binary_positive_reference` or set
  `reference_sign = -1`.

**returning structure**: A list containing:

- `scRNA_data`: A Seurat object after screening, with the group labels
  stored in the `PhenoMapR` column (`"Positive"`, `"Negative"`,
  `"Neutral"` or `NA`) and the per-cell scores in `PhenoMapR_*` columns
  of `meta.data`. The parameters used are recorded in
  `scRNA_data@misc$PhenoMapR`.
- `phenomap_scores`: Per-cell scores returned by
  `PhenoMapR::PhenoMap()`.
- `reference`: The gene-level reference z-scores derived from the bulk
  data.

#### 2.F Merge screening results

If you have performed multiple screening methods on **the same**
single-cell data, you can use the `MergeResult` to merge the screening
results of these methods. The Seurat object or a results list from
`Screen` is accepted.

``` r

merged_seurat <- MergeResult(
  your_scissor_result,
  your_scPAS_result,
  your_scAB_result,
  your_scPP_result,
  your_DEGAS_result
  # # * Add more if you want
  # ,your_LP_SGL_result,
  # your_PIPET_result,
  # your_SIDISH_result,
  # your_SCIPAC_result,
  # your_DEGASv2_result,
  # your_TiRank_result,
  # your_PhenoMapR_result
)

# * mixed input form is alse supported

merged_seurat <- MergeResult(
  your_scissor_result$scRNA_data,
  your_scPAS_result$scRNA_data,
  your_scAB_result,
  your_scPP_result,
  your_DEGAS_result
  # # * Add more if you want
  # ,your_LP_SGL_result$scRNA_data,
  # your_PIPET_result
)
```

This function merges all slots from the input Seurat objects. Note that
intermediate data (e.g., `scissor_result$reliability.test` or
`scab_result$scAB_result`) **will not be retained** during the merge.
While the function can technically combine different single-cell
datasets, it is specifically designed for merging replicates or batches
of the **same** single-cell RNA-seq experiment; using it for
heterogeneous data may lead to subtle (and potentially hard-to-detect)
errors (You may consider `SCIntegrate` for this purpose).

**returning structure**:

A Seurat object with all merged slots.

------------------------------------------------------------------------

### 3. Visualization

Here we provide some visualization methods for the screening results.
Considering that many people have different needs for data
visualization, `SigBridgeR` hardly provides visualization (except for
fraction plot and upset plot, because we provide some statistic results
for them).

If you have any suggestions for visualization, please submit an issue or
pull request.

#### 3.1 UMAP

**example**:

Suppose you have performed all algorithm screening on your Seurat object
and wish to examine the distribution across different celltypes and
patient, you may reference and use the following code:

``` r

# * Glance your seurat object
ggplot2::autoplot(merged_seurat, group.by = "scissor")
```

``` r

library(zeallot)
# library(Seurat)
# library(patchwork)
# library(purrr)

c(
  celltype_umap,
  patient_umap,
  scissor_umap,
  scab_umap,
  scpas_umap,
  scpp_umap,
  degas_umap,
  # * Add more if you want
) %<-%
  purrr::map(
    # make sure these column names exist
    c("celltype", "patient", "scissor", "scAB", "scPAS", "scPP", "DEGAS"),
    ~ Seurat::DimPlot(
      merged_seurat,
      group.by = .x,
      pt.size = 0.05,
      reduction = "umap"
    ) +
      ggplot2::ggtitle(.x)
  )

# * Show
umaps <- celltype_umap +
  patient_umap +
  scissor_umap +
  scab_umap +
  scpas_umap +
  scpp_umap +
  degas_umap +
  # * Add more if you want
  patchwork::plot_layout(ncol = 3)

umaps
```

This will generate seven UMAP plots separately.

Or suppose you have performed `scPAS` screening on your Seurat object
and want to visualize the distribution of prediction confidence scores,
you may reference and use the following code:

``` r

library(zeallot)
library(patchwork)
# library(Seurat)
# library(purrr)

c(scPAS_Pvalue_umap, scPAS_NRS_umap) %<-%
  purrr::map(
    c("scPAS_Pvalue", "scPAS_NRS"),
    ~ Seurat::FeaturePlot(
      object = merged_seurat,
      features = .x,
    ) +
      ggplot2::ggtitle(.x) +
      ggplot2::theme(legend.position = "right")
  )

# * Show
scPAS_Pvalue_umap | scPAS_NRS_umap
```

This will generate two plots, one for each feature specified in
`feature`.

**helpful documentation**:

[Seurat::DimPlot -
https://satijalab.org/seurat/reference/DimPlot.html](https://satijalab.org/seurat/reference/DimPlot.html)

[Seurat::FeaturePlot -
https://satijalab.org/seurat/reference/FeaturePlot.html](https://satijalab.org/seurat/reference/FeaturePlot.html)

#### 3.2 Stack bar plot

Key parameters for `ScreenFractionPlot`:

- `seurat_obj`: A Seurat object after screening.
- `group_by`: Used to specify the column of the meta.data in
  `seurat_obj`. The plot results will be grouped by this parameter.
- `screen_type`: Screening algorithm used before. (case-sensitive, e.g.,
  “scissor” for Scissor results)
- `show_null`: Logical, whether to show groups with zero cells (default:
  FALSE).
- `plot_color` Custom color palette (named vector format):
  - Required names: “Positive”, “Negative”, “Neutral”, “Other”
  - Default: c(“Neutral”=“#CECECE”, “Other”=“#CECECE”,
    “Positive”=“#ff3333”, “Negative”=“#386c9b”)

Suppose you have already performed the `scPAS` algorithm screening on
your Seurat object, and you want to check the proportion of positive
cells across different patients. You can refer to and use the following
code:

``` r

plot <- ScreenFractionPlot(
  screened_seurat = merged_seurat,
  group_by = "seurat_clusters", # grouping basis for the x-axis
  screen_type = "scPAS",
  plot_title = "scPAS Screening Results"
)
plot
```

If you have performed multiple screening methods and already merged the
results, you can use the following code:

``` r

plot <- ScreenFractionPlot(
  screened_seurat = merged_seurat,
  group_by = "seurat_clusters", # grouping basis for the x-axis
  screen_type = c("scPAS", "scissor"), # multiple screening results
  plot_title = "Screening Results" #  A screen_type prefix will be added to the current plot title
)
plot$combined_plot
```

The order of the groups is determined by the proportion of **Positive**
cells within each group.

**returning structure**:

If a single screen_type is specified

- `stats`: A data frame containing the proportion of positive cells for
  each group.
- `plot`: A ggplot2 object.

If multiple screen_types are specified

- `stats`: A list containing data frames containing the proportion of
  positive cells for each group.
- `plot`: A list containing each ggplot2 objects.
- `combined_plot`: A ggplot2 object containing all the plots (2\*2
  grid).

#### 3.3 Venn diagram

`ggVennDiagram` is used to generate a Venn diagram for the screening
results. Suppose you have performed some of the screening algorithms on
your Seurat object, and you want to check the overlap of the cells
selected by each algorithm. You can refer to and use the following code:

**example**:

``` r

library(ggVennDiagram)

# # * If you have merged the results, you can use the following code instead:
# c(scissor_pos, scab_pos, scpas_pos, scpp_pos, degas_pos) %<-%
#     purrr::map(
#         c("scissor", "scAB", "scPAS", "scPP", "DEGAS"),
#         ~ colnames(merged_seurat)[
#             which(merged_seurat[[.x]] == "Positive")
#         ]
#     )

# * get the cell vectors
scissor_pos <- colnames(scissor_result$scRNA_data)[
  which(scissor_result$scRNA_data$scissor == "Positive")
]
scab_pos <- colnames(scab_result$scRNA_data)[
  which(scab_result$scRNA_data$scAB == "Positive")
]
scpas_pos <- colnames(scpas_result$scRNA_data)[
  which(scpas_result$scRNA_data$scPAS == "Positive")
]
scpp_pos <- colnames(scissor_result$scRNA_data)[
  which(scpp_result$scRNA_data$scPP == "Positive")
]
degas_pos <- colnames(degas_result$scRNA_data)[
  which(degas_result$scRNA_data$DEGAS == "Positive")
]


all_cells <- colnames(your_seurat_obj)

# * create a list of cell vectors
pos_venn <- list(
  scissor = scissor_pos,
  scpas = scpas_pos,
  scab = scab_pos,
  scpp = scpp_pos,
  degas = degas_pos,
  all_cells = all_cells
  # * you can add more groups here
)

set.seed(123)

venn_plot <- ggVennDiagram::ggVennDiagram(
  x = pos_venn,
  # * the labels of each group to be shown on the diagram
  category.names = c(
    "Scissor",
    "scPAS",
    "scAB",
    "scPP",
    "DEGAS",
    "All cells"
  ),
  # * the colors of each group
  set_color = c(
    "#a33333ff",
    "#37ae00ff",
    "#2a2a94ff",
    "#9c8200ff",
    "#bb14adff",
    "#008383ff"
  )
) +
  ggplot2::scale_fill_gradient(low = "white", high = "#ffb6b6ff") +
  ggplot2::ggtitle("Screening Venn Diagram")
```

**helpful link**:

[ggVennDiagram -
https://gaospecial.github.io/ggVennDiagram/](https://gaospecial.github.io/ggVennDiagram/)

#### 3.4 Upset plot

If too many screening meyhods are selected, the number of intersections
among cells screened by different methods will also increase. In this
case, using an upset plot is more intuitive and neat than a Venn
diagram.

`ggupset` is used for visualizing upset plot.

Key parameters for `ScreenUpset`:

- `screened_seurat`: A Seurat object after screening.
- `screen_type`: Screening algorithm used before. (case-sensitive, e.g.,
  “scissor” for Scissor results)
- `n_intersections` : Number of intersections to display in the plot.
  Default: 20.
- `x_lab` : Label for the x-axis. Default: “Screen Set Intersections”.
- `y_lab` : Label for the y-axis. Default: “Number of Cells”.
- `title` : Plot title. Default: “Cell Counts Across Screen Set
  Intersections”.
- `bar_color` : Color for the bars in the plot. Default: “#4E79A7”.
- `combmatrix_point_color` : Color for points in the combination matrix.
  Default: “black”.
- `...` : Additional arguments passed to
  [`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html)
  for customizing the plot appearance.

``` r

upset <- ScreenUpset(screened_seurat = merged_seurat)
upset$plot
```

**returning structure**:

- `plot`: A ggplot2 object.
- `stats`: A data frame containing the numbers of positive cells for
  each group.

**helpful link**:

[ggupset -
https://github.com/const-ae/ggupset](https://github.com/const-ae/ggupset)

------------------------------------------------------------------------

### 4. Example (last updated: 2025-7-31)

#### 4.1 Survival-associated cell screening

Here we use the example data LUAD to demonstrate how to use the
functions in `SigBridgeR` to screen cells associated with phenotype.

``` r

library(SigBridgeR)
library(zeallot)

# * load the example data
c(seurat, bulk, pheno) %<-% LoadRefData(data_type = "survival")

seurat
# An object of class Seurat
# 15924 features across 1077 samples within 1 assay
# Active assay: RNA (15924 features, 2000 variable features)
#  3 layers present: counts, data, scale.data
#  2 dimensional reductions calculated: pca, umap

dim(bulk)
# [1] 4071  506
bulk[1:6, 1:6] # already log2 transformed
#         TCGA-69-7978 TCGA-62-8399 TCGA-78-7539 TCGA-73-4658 TCGA-44-6775 TCGA-44-2655
# HIF3A         4.2598      11.6239       9.1362       5.0288       4.0573       5.5335
# RTN4RL2       8.2023       5.5819       3.5365       7.4156       7.7107       5.3257
# HMGCLL1       2.7476       5.8513       3.8334       3.6447       2.9188       4.8820
# LRRTM1        0.0000       0.4628       4.7506       6.8005       7.7819       2.2882
# GRIN1         6.6074       5.4257       4.9563       7.3510       3.5361       3.3311
# LRRTM3        1.7458       2.0092       0.0000       1.4468       0.0000       0.0000
nrow(pheno)
# [1] 506
head(pheno)
#               time status
# TCGA-69-7978  4.40      0
# TCGA-62-8399 88.57      0
# TCGA-78-7539 25.99      0
# TCGA-73-4658 52.56      1
# TCGA-44-6775 23.16      0
# TCGA-44-2655 43.50      0
```

This single-cell RNA data is from humans. It has already been
preprocessed into a **Seurat** object (`seurat`) with QC filtering,
normalization, clustering, and UMAP included, so it can be passed
directly to
[`Screen()`](https://wanglabcsu.github.io/sigbridger/reference/Screen.md).

First, we use `scissor` to screen cells associated with survival.

``` r

scissor_result <- Screen(
  matched_bulk = bulk,
  sc_data = seurat,
  phenotype = pheno,
  label_type = "survival",
  phenotype_class = "survival",
  screen_method = "Scissor",
  alpha = 0.05
)

# i [2025/09/08 17:03:20] Scissor start...
# i [2025/09/08 17:03:20] Start from raw data...
# i Using "RNA_snn" graph for network.
# i [2025/09/08 17:03:20] Normalizing quantiles of data...
# i [2025/09/08 17:03:20] Subsetting data...
# i [2025/09/08 17:03:21] Calculating correlation...
# ----------------------------------------------------------------------------------------------------
# Five-number summary of correlations:

#         0%        25%        50%        75%       100%
# -0.2342323  0.0594385  0.1118708  0.1642065  0.5250605
# ----------------------------------------------------------------------------------------------------
# i [2025/09/08 17:03:21] Perform cox regression on the given clinical outcomes:
# v Statistics data saved to `Scissor_inputs.RData`.
# i [2025/09/08 17:03:22] Screening...
#
# -- At alpha = 0.05 --
#
# Scissor identified 265 Scissor+ cells and 73 Scissor- cells.
# The percentage of selected cell is: 30.924%

table(scissor_result$scRNA_data$scissor)
# Negative  Neutral Positive
#       73      755      265
```

You will see an additional “Scissor_inputs.RData” in the working
directory. This is the intermediate data generated by the Scissor
method, which we can use to save running time ( You can also set
`path2load_scissor_cache=NULL` to suppress the saving of intermediate
files ). Meanwhile, we set `reliability_test=TRUE`, which will run an
additional reliability test.

``` r

scissor_result <- Screen(
  matched_bulk = bulk, # doesn't need to be provided Since the intermediate data is provided
  sc_data = seurat,
  phenotype = pheno_ok, # doesn't need to be provided Since the intermediate data is provided
  label_type = "survival",
  phenotype_class = "survival",
  screen_method = "Scissor",
  alpha = 0.05,
  path2load_scissor_cache = "Scissor_inputs.RData",
  reliability_test = TRUE
)
# i [2025/09/08 16:07:48] Scissor start...
# i [2025/09/08 16:07:48] Loading data from `Scissor_inputs.RData`...
# i [2025/09/08 16:07:48] Screening...
# [1] "alpha = 0.05"
# [1] "Scissor identified 249 Scissor+ cells and 245 Scissor- cells."
# [1] "The percentage of selected cell is: 45.197%"

# --------------------------------------------------------------------------------
# i [2025/09/08 16:07:54] Start reliability test

# Attaching package: 'survival'

# The following object is masked from 'package:future':

#     cluster

# [1] "|**************************************************|"
# [1] "Perform cross-validation on X with true label"
# Finished!
# [1] "|**************************************************|"
# [1] "Perform cross-validation on X with permutated label"
# Finished!
# [1] "Test statistic = 0.590"
# [1] "Reliability significance test p = 0.000"
# v [2025/09/08 16:10:54] reliability test: Done

scissor_result$reliability_result$statistic
# [1] 0.5899587

scissor_result$reliability_result$p
# [1] 0

scissor_result$reliability_result$c_index_test_real
#  [1] 0.6038544 0.5022222 0.6050725 0.6279391 0.5064935 0.6033520 0.6769231 0.6453089 0.4968421 0.6315789

scissor_result$reliability_result$c_index_test_back |>
  unlist() |>
  matrix(nrow = 10)
#            [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]      [,8]      [,9]     [,10]
#  [1,] 0.5594714 0.4943820 0.5379747 0.5583658 0.5527728 0.6146435 0.5130597 0.5701275 0.4691943 0.5177305
#  [2,] 0.5735608 0.4726891 0.4146341 0.5020661 0.5354331 0.6840149 0.5540275 0.4990758 0.5324484 0.5497382
#  [3,] 0.5960145 0.5157116 0.6191446 0.5091912 0.5659656 0.5870370 0.5090580 0.4247788 0.5734266 0.5539715
#  [4,] 0.6210191 0.6563147 0.5581948 0.4675615 0.4688222 0.5308219 0.5100402 0.6709091 0.6155779 0.5607143
#  [5,] 0.4623656 0.5716695 0.4920441 0.5403727 0.5555556 0.5911215 0.5738499 0.6189258 0.5871560 0.6084337
#  [6,] 0.5868794 0.5127660 0.7416880 0.5366876 0.5417440 0.5633803 0.5161943 0.4830372 0.4732965 0.5740132
#  [7,] 0.4975610 0.5751503 0.4801902 0.5588972 0.4940898 0.5723370 0.5788382 0.6171875 0.5884413 0.5445545
#  [8,] 0.6331361 0.5970516 0.5473888 0.6274131 0.5159705 0.5000000 0.5299760 0.4905660 0.5277778 0.5027322
#  [9,] 0.6390805 0.6548913 0.5049310 0.6299639 0.6385135 0.5964392 0.5381605 0.7462687 0.5185185 0.6585859
# [10,] 0.5207373 0.5000000 0.6051780 0.4932127 0.6892430 0.4786517 0.5619266 0.6614583 0.6502242 0.5333333
```

Next, we use scPAS, scAB, scPP, LP_SGL and PIPET to do the same
screening. Generally you only need to change the `screen_method`, as
long as you have not specified any particular parameters.

``` r

scpas_result <- Screen(
  matched_bulk = bulk,
  sc_data = seurat,
  phenotype = pheno,
  label_type = "survival",
  phenotype_class = "survival",
  screen_method = "scPAS",
  alpha = 0.05
)
# i [2025/10/20 16:24:27] Start scPAS screening.
# i [2025/10/20 16:24:29] Quantile normalization of bulk data.
# i [2025/10/20 16:24:29] Extracting single-cell expression profiles...
# i [2025/10/20 16:24:29] Constructing a gene-gene similarity by single cell data...
# Building SNN based on a provided distance matrix
# Computing SNN
# i [2025/10/20 16:24:30] Optimizing the network-regularized sparse regression model...
# i [2025/10/20 16:24:30] Perform cox regression on the given phenotypes...
#
# -- At alpha = 0.05 --
#
# lambda = 0.776168989003421
# scPAS identified 59 risk+ features and 61 risk- features.
# The percentage of selected feature is: 13.73%
# i [2025/10/20 16:24:49] Calculating quantified risk scores...
# i [2025/10/20 16:24:49] Qualitative identification by permutation test program with 2000 times random perturbations...
# v [2025/10/20 16:24:50] scPAS screening done.

table(scpas_result$scRNA_data$scPAS)

# Negative  Neutral Positive
#       5     1085        3
```

As you can see, due to differences in data and algorithms, not every
screening algorithm is able to screen out cells. You can adjust the
corresponding parameters, e.g. change the `alpha` to `NULL`, this will
make scPAS iterate alpha until the result is significant (or judged as
having no significant cell subpopulations). See also [3.2 (Option B)
scPAS Screening](#id_32-option-b-scpas-screening)

``` r

scpas_result <- Screen(
  matched_bulk = bulk,
  sc_data = seurat,
  phenotype = pheno,
  label_type = "survival",
  phenotype_class = "survival",
  screen_method = "scPAS",
  alpha = NULL,
  cutoff = 0.2
)
# i [2025/10/20 16:43:31] Start scPAS screening.
# i [2025/10/20 16:43:32] Quantile normalization of bulk data.
# i [2025/10/20 16:43:32] Extracting single-cell expression profiles...
# i [2025/10/20 16:43:32] Constructing a gene-gene similarity by single cell data...
# Building SNN based on a provided distance matrix
# Computing SNN
# i [2025/10/20 16:43:33] Optimizing the network-regularized sparse regression model...
# i [2025/10/20 16:43:33] Perform cox regression on the given phenotypes...
#
# -- At alpha = 0.001 --
#
# lambda = 7.61825638357188
# scPAS identified 315 risk+ features and 366 risk- features.
# The percentage of selected feature is: 77.918%
#
# -- At alpha = 0.005 --
#
# lambda = 3.68743152046651
# scPAS identified 193 risk+ features and 231 risk- features.
# The percentage of selected feature is: 48.513%
#
# -- At alpha = 0.01 --
#
# lambda = 2.55333682907113
# scPAS identified 137 risk+ features and 158 risk- features.
# The percentage of selected feature is: 33.753%
#
# -- At alpha = 0.05 --
#
# lambda = 0.776168989003421
# scPAS identified 59 risk+ features and 61 risk- features.
# The percentage of selected feature is: 13.73%
# i [2025/10/20 16:44:55] Calculating quantified risk scores...
# i [2025/10/20 16:44:55] Qualitative identification by permutation test program with 2000 times random perturbations...
# v [2025/10/20 16:44:57] scPAS screening done.
```

Generally speaking, a larger alpha means looser screening. However, it
still won’t result in the occurence of too many false-positive cells.
You can also adjust the `cutoff` parameter as you like.

``` r

table(scpas_result$scRNA_data$scPAS)
# Negative  Neutral Positive
#        5     1085        3
```

Now we use scAB, scPP, DEGAS, LP_SGL and PIPET to screen cells.

``` r

scab_result <- Screen(
  matched_bulk = bulk,
  sc_data = seurat,
  phenotype = pheno,
  label_type = "survival",
  phenotype_class = "survival",
  screen_method = "scAB"
)
# i [2025/09/08 17:04:51] Start scAB screening.
# i  Using "RNA_snn" graph for network.
# i [2025/09/08 17:04:52] Selecting K...
# i [2025/09/08 17:06:15] Run NMF with phenotype and cell-cell similarity regularization at K = 3.
# i [2025/09/08 17:06:19] Screening cells...
# i [2025/09/08 17:06:19] scAB screening done.

table(scab_result$scRNA_data$scAB)
#    Other Positive
#     1018       75
```

``` r

scpp_result <- Screen(
  matched_bulk = bulk,
  sc_data = seurat,
  phenotype = pheno,
  label_type = "survival",
  phenotype_class = "survival",
  screen_method = "scPP"
)
# i [2025/09/08 17:00:28] Start scPP screening.
# i [2025/09/08 17:00:28] Finding markers...
# Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights,  :
#   Loglik converged before variable  1 ; coefficient may be infinite.
# i [2025/09/08 17:00:52] Screening...
# Genes in the gene sets NOT available in the dataset:
#   gene_pos:   13 (6% of 230)
#   gene_neg:   54 (12% of 446)
# There are no genes significantly upregulated in Phenotype- compared to Phenotype+.
# v [2025/09/08 17:00:54] scPP screening done.

table(scpp_result$scRNA_data$scPP)
# Negative  Neutral Positive
#       52      993       48
```

``` r

# I recommend running this code in the background.
degas_result <- Screen(
  matched_bulk = bulk,
  sc_data = seurat,
  phenotype = pheno,
  label_type = "This_is_a_DEGAS_test",
  phenotype_class = "survival",
  screen_method = "DEGAS"
)
# i [2025/10/09 18:41:00] Starting DEGAS Screen
# i [2025/10/09 18:41:01] Setting up Environment...
# i [2025/10/09 18:41:08] Training DEGAS model...
# i [2025/10/09 18:41:08] 3-layer DenseNet BlankCox DEGAS model
# i [2025/10/09 18:41:10] Python check passed, using Python 3.9.15
# i [2025/10/09 18:41:10] Training...
###########################
# Many output from python #
###########################
# i [2025/10/10 17:35:04] Predicting and Labeling...
# i [2025/10/10 17:35:04] Labeling screened cells...
# i [2025/10/10 17:35:04] Searching for survival-associated cells...
# i Scores over 0.499 are considered `Positive`.
# i [2025/10/10 17:35:04] DEGAS Screen done.

table(degas_result$scRNA_data$DEGAS)
#    Other Positive
#     1038       55
```

``` r

# I recommend running this code in the background.
lpsgl_result <- Screen(
  matched_bulk = bulk,
  sc_data = seurat,
  phenotype = pheno,
  label_type = "LP_SGL",
  phenotype_class = "survival",
  screen_method = "LP_SGL"
)
# i [2025/11/22 21:14:29] Starting LP-SGL screening analysis
# i [2025/11/22 21:14:29] Fetch graph from Seurat object
# i [2025/11/22 21:14:29] Run Leiden clustering
# i [2025/11/22 21:14:30] Calculating correlation matrix...
# i [2025/11/22 21:14:33] Fitting SGL model with alpha = 0.5, this may take a while
# i [2025/11/22 21:14:49] Running 5-fold cross-validation
# i [2025/11/22 21:17:31] Optimal lambda index: 20 (error = 1380.7537932881)
# i [2025/11/22 21:17:31] LP-SGL screening completed

table(lpsgl_result$scRNA_data$LP_SGL)
# Negative  Neutral Positive
#       77      778      233
```

Note that PIPET does not support survival phenotype. We will show how to
use PIPET in [Section
5.3](#id_53-binarized-phenotype-associated-cell-screening)

After these algorithms have been run, the four sets of data can be
merged since screening methods performed on the same data.

``` r

screen_result <- MergeResult(
  scissor_result,
  scpas_result,
  scab_result,
  scpp_result,
  degas_result,
  lpsgl_result
)
# v Successfully merged 6 objects.

class(screen_result)
# [1] "Seurat"
# attr(,"package")
# [1] "SeuratObject"

colnames(screen_result[[]])
#  [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "test_col"        "percent.rp"      "percent.mt"      "RNA_snn_res.0.1" "seurat_clusters"
#  [9] "scissor"         "scAB"            "scAB_Subset1"    "Subset1_loading" "scAB_Subset2"    "Subset2_loading" "scPAS_RS"        "scPAS_NRS"
# [17] "scPAS_Pvalue"    "scPAS_FDR"       "scPAS"           "scPP_AUCup"      "scPP_AUCdown"    "scPP"            "DEGAS"           "LP_SGL"
```

Finally, we can visualize the screening results. Let’s start with a Venn
diagram to see the situation.

``` r

library(ggVennDiagram)
library(zeallot)

# * color palette
set.seed(123)
my_colors <- RColorBrewer::brewer.pal(
  length(unique(screen_result$seurat_clusters)),
  "Set1"
)

c(scissor_pos, scab_pos, scpas_pos, scpp_pos, degas_pos, lpsgl_pos) %<-%
  purrr::map(
    c("scissor", "scAB", "scPAS", "scPP", "DEGAS", "LP_SGL"),
    ~ colnames(screen_result)[
      which(screen_result[[.x]] == "Positive")
    ]
  )

all_cells <- colnames(screen_result)

# * create a list of cell vectors
pos_venn <- list(
  scissor = scissor_pos,
  scpas = scpas_pos,
  scab = scab_pos,
  scpp = scpp_pos,
  degas = degas_pos,
  lp_sgl = lpsgl_pos,
  all_cells = all_cells
)


venn <- ggVennDiagram::ggVennDiagram(
  x = pos_venn,
  # * the labels of each group to be shown on the diagram
  category.names = c(
    "Scissor",
    "scPAS",
    "scAB",
    "scPP",
    "DEGAS",
    "LP_SGL",
    "All cells"
  ),
  # * the colors of each group
  set_color = c(
    "#a33333ff",
    "#37ae00ff",
    "#0000f5ff",
    "#d4b100ff",
    "#e600eeff",
    "#59108aff",
    "#008383ff"
  ),
  label_geom = "text"
) +
  ggplot2::scale_fill_gradient(low = "white", high = "#ffb6b6ff") +
  ggplot2::ggtitle("Screening Venn Diagram")

venn
```

[![venn](articles/example_figures/venn.png)](https://wanglabcsu.github.io/sigbridger/articles/vignettes/example_figures/venn.png)

When there are too many sets, the visualization effect of the Venn
diagram is not very good. We can use a set plot instead. Only positive
cells are shown in the set plot.

``` r

upset <- ScreenUpset(
  screened_seurat = screen_result,
  screen_type = c("scissor", "scPAS", "scAB", "scPP", "DEGAS", "LP_SGL"),
  n_intersections = 40
)

# * show the cell numbers of each set
head(upset$stats)
# # A tibble: 6 x 3
#   intersection    sets         count
#   <chr>           <named list> <dbl>
# 1 scissor         <chr [1]>      265
# 2 scPAS           <chr [1]>        2
# 3 scAB            <chr [1]>       75
# 4 scPP            <chr [1]>       49
# 5 DEGAS           <chr [1]>       55
# 6 LP_SGL          <chr [1]>      233

# ggplot2::ggsave(
#   "vignettes/vignettes/example_figures/upset.png",
#   plot = upset$plot,
#   width = 10,
#   height = 10
# )
```

[![upset](articles/example_figures/upset.png)](https://wanglabcsu.github.io/sigbridger/articles/vignettes/example_figures/upset.png)

A bar chart showing proportions can also be used to examine the
screening results. Since the example data does not have sample metadata,
we have created a fictional `Sample` column.

``` r

set.seed(123)
# * fictional sample column
screen_result$Sample <- sample(
  paste0("Sample", 1:10),
  ncol(screen_result),
  replace = TRUE
)

table(screen_result$Sample)
#  Sample1 Sample10  Sample2  Sample3  Sample4  Sample5  Sample6  Sample7  Sample8  Sample9
#       98      117       96      119       95      100      101      131      115      121

fraction_list <- ScreenFractionPlot(
  screened_seurat = screen_result,
  group_by = "Sample",
  screen_type = c("scissor", "scPP", "scAB", "scPAS", "DEGAS", "LP_SGL"),
  show_null = FALSE,
  plot_color = NULL,
  show_plot = TRUE
)
# Creating plots for 6 screen types...

# ggplot2::ggsave(
#   "vignettes/vignettes/example_figures/fraction.png",
#   plot = fraction_list$combined_plot,
#   width = 10,
#   height = 10
# )
```

[![fraction](articles/example_figures/fraction.png)](https://wanglabcsu.github.io/sigbridger/articles/vignettes/example_figures/fraction.png)

As you see, the `label_type`s set in function `Screen` are shown in the
legend of each plot.

The `fraction_list` contains the statistical data and charts for each
screening algorithm. The title of each plot will be appended with the
method name as an identifier.

fraction_list  
\|– stats  
\| \|– scissor  
\| \|– scAB  
\| \|– scPAS  
\| \|– scPP  
\| \|– DEGAS  
\| -- LP_SGL  
\|– plots  
\| \|– scissor  
\| \|– scAB  
\| \|– scPAS  
\| \|– scPP  
\| \|– DEGAS  
\| -- LP_SGL  
-- combined_plot \# show 6 plots in one plot

UMAP is the most commonly used type of plot in academic literature.

``` r

library(patchwork)
library(zeallot)

my_palette <- RColorBrewer::brewer.pal(
  n = length(unique(screen_result$Sample)),
  "Set1"
)

sample_umap <- Seurat::DimPlot(
  screen_result,
  group.by = "Sample",
  pt.size = 1.2,
  alpha = 0.8,
  reduction = "umap",
  cols = my_palette
) +
  ggplot2::ggtitle("Sample")


c(
  scissor_umap,
  scab_umap,
  scpas_umap,
  scpp_umap,
  degas_umap,
  lpsgl_umap
) %<-%
  purrr::map(
    c("scissor", "scAB", "scPAS", "scPP", "DEGAS", "LP_SGL"), # make sure these column names exist
    ~ Seurat::DimPlot(
      screen_result,
      group.by = .x,
      pt.size = 1.2,
      alpha = 0.8,
      reduction = "umap",
      cols = c(
        "Neutral" = "#CECECE",
        "Other" = "#CECECE",
        "Positive" = "#ff3333",
        "Negative" = "#386c9b"
      )
    ) +
      ggplot2::ggtitle(.x)
  )

# * Show
umaps <- sample_umap +
  scissor_umap +
  scab_umap +
  scpas_umap +
  scpp_umap +
  degas_umap +
  lpsgl_umap +
  patchwork::plot_layout(ncol = 2)

umaps

# ggplot2::ggsave(
#      "vignettes/vignettes/example_figures/umaps.png",
#      plot = umaps,
#      width = 10,
#      height = 10
# )
```

[![umaps](articles/example_figures/umaps.png)](https://wanglabcsu.github.io/sigbridger/articles/vignettes/example_figures/umaps.png)

All the analysis has been completed, and we can save the data for future
use.

``` r

SeuratObject::SaveSeuratRds(
  object = screen_result,
  file = "screened_result.rds"
)
```

or in `.h5ad` format:

``` r

anndataR::write_h5ad(
  object = screen_result,
  path = "screened_result.h5ad",
  compression = "gzip"
)
```

#### 4.2 Continuous phenotype associated cell screening

Generally speaking, the process is the same as described in [5.1
Survival-associated cell
screening](#id_51-survival-associated-cell-screening). Here, only the
preprocessing of continuous phenotypic data is introduced.

``` r

library(SigBridgeR)
library(zeallot)
library(Seurat)
setwd(here::here())

# * load the example data
c(seurat, bulk, pheno) %<-% LoadRefData(data_type = "continuous")

seurat
# An object of class Seurat
# 15911 features across 1077 samples within 1 assay
# Active assay: RNA (15911 features, 2000 variable features)
#  3 layers present: counts, data, scale.data
#  2 dimensional reductions calculated: pca, umap
dim(bulk)
# [1] 4106  289
bulk[1:6, 1:6]
#        TCGA-AZ-6599-01 TCGA-AA-3655-01 TCGA-A6-6137-01 TCGA-CK-4952-01 TCGA-A6-5657-01 TCGA-AD-6963-01
# HIF3A           2.3437          2.0858          6.0759          1.9506          5.4777          4.4634
# CAMK4           4.9331          2.3709          4.1387          1.1557          4.1746          3.2363
# RNF112          2.4817          2.4947          3.5941          2.3486          4.9185          1.4621
# SPN             5.6704          6.8577          8.0598          5.0049          7.6076          7.3960
# LRRTM1          1.6031          0.9465          1.9142          0.0000          3.2523          0.0000
# GRIN1           6.4944          4.3225          2.8073          7.3460          4.5000          3.1816

# * A named vector
head(pheno)
# TCGA-AZ-6599-01 TCGA-AA-3655-01 TCGA-A6-6137-01 TCGA-CK-4952-01 TCGA-A6-5657-01 TCGA-AD-6963-01
#             178              65              91             206              63              67
```

If your phenotype is a `data.frame`, try this to convert it to a
`named vector`:

``` r

pheno <- setNames(your_data.frame$continuous_numeric, your_data.frame$sample)
```

#### 4.3 Binarized phenotype associated cell screening

This process is also the same as described in [5.1 Survival-associated
cell screening](#id_51-survival-associated-cell-screening). Here, only
the preprocessing of binary phenotypic data is introduced.

``` r

library(SigBridgeR)
library(zeallot)
library(Seurat)

# * load the example data
c(seurat, bulk, pheno) %<-% LoadRefData(data_type = "binary")

seurat
# An object of class Seurat
# 15911 features across 1077 samples within 1 assay
# Active assay: RNA (15911 features, 2000 variable features)
#  3 layers present: counts, data, scale.data
#  2 dimensional reductions calculated: pca, umap
dim(bulk)
# [1] 4106  434
bulk[1:6, 1:6]
#        TCGA-CA-5256-01 TCGA-AZ-6599-01 TCGA-AA-3655-01 TCGA-A6-6137-01 TCGA-CK-4952-01 TCGA-A6-5657-01
# HIF3A           3.7172          2.3437          2.0858          6.0759          1.9506          5.4777
# CAMK4           3.0698          4.9331          2.3709          4.1387          1.1557          4.1746
# RNF112          1.3702          2.4817          2.4947          3.5941          2.3486          4.9185
# SPN             5.5207          5.6704          6.8577          8.0598          5.0049          7.6076
# LRRTM1          3.2408          1.6031          0.9465          1.9142          0.0000          3.2523
# GRIN1           3.0698          6.4944          4.3225          2.8073          7.3460          4.5000

# * A named vector
head(pheno)
# TCGA-CA-5256-01 TCGA-AZ-6599-01 TCGA-AA-3655-01 TCGA-A6-6137-01 TCGA-CK-4952-01 TCGA-A6-5657-01
#               1               1               1               1               1               1
```

If your phenotype is a `data.frame`, your binary variable is stored in
the “data” column, categorized as ‘Tumor’ and ‘Normal’, we will assign
‘Tumor’ a value of 1 and ‘Normal’ a value of 0. In this way, cells
screened as **Positive** will be associated with ‘Tumor’. Try this to
convert it to a `named vector`:

``` r

pheno <- mutate(
  pheno,
  data = dplyr::case_when(
    data == "Tumor" ~ 1,
    data == "Normal" ~ 0
  )
)
pheno <- setNames(pheno$data, pheno$Sample)
```

**A PIPET use case**

Previously, it was noted that PIPET cannot use patient survival data as
the phenotype; here, we demonstrate using a binary phenotype instead.

``` r

pipet_result <- Screen(
  matched_bulk = bulk,
  sc_data = seurat,
  phenotype = pheno,
  phenotype_class = "binary",
  screen_method = "PIPET",
  label_type = "PIPET"
)
#i [2025/12/26 06:18:28] Starting PIPET screen
#i [2025/12/26 06:18:28] Creating marker genes from bulk data...
#i [2025/12/26 06:18:28] Fitting model with limma
#v [2025/12/26 06:18:28] Created 2795 marker genes
#i [2025/12/26 06:18:28] Running PIPET correlation analysis...
#i [2025/12/26 06:18:28] The classification of markers is: group_0: 945 and group_1: 477
#i [2025/12/26 06:18:28] Normalize count data with CPM and log1p
#i [2025/12/26 06:18:29] Scale features with z-score normalization
#i [2025/12/26 06:19:49] Organize the computed results
#v [2025/12/26 06:19:49] PIPET screening done.

table(pipet_result$scRNA_data$PIPET)
# Negative  Neutral Positive
#     165      708      220
```

------------------------------------------------------------------------

### 5. Troubleshooting

View
Troubleshooting[(https://github.com/WangX-Lab/SigBridgeR/wiki/Troubleshooting)](https://wanglabcsu.github.io/SigBridgeR/articles/Troubleshooting.html)

Session information:

    ## R version 4.6.2 (2026-10-29)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 26.04.1 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] SigBridgeR_4.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] pbapply_1.7-5         gridExtra_2.3.1       rlang_1.3.0          
    ##   [4] magrittr_2.0.5        RcppAnnoy_0.0.23      otel_0.3.0           
    ##   [7] matrixStats_1.5.0     ggridges_0.5.7        compiler_4.6.2       
    ##  [10] png_0.1-9             systemfonts_1.3.4     vctrs_0.7.3          
    ##  [13] reshape2_1.4.5        stringr_1.6.0         pkgconfig_2.0.3      
    ##  [16] fastmap_1.2.0         XVector_0.52.0        rmarkdown_2.32       
    ##  [19] ragg_1.5.2            purrr_1.2.2           xfun_0.62            
    ##  [22] cachem_1.1.0          beachmat_2.28.0       jsonlite_2.0.0       
    ##  [25] DelayedArray_0.38.2   chk_0.11.0            irlba_2.4.1          
    ##  [28] parallel_4.6.2        cluster_2.1.8.3       R6_2.6.1             
    ##  [31] ica_1.0-3             bslib_0.12.0          stringi_1.8.9        
    ##  [34] RColorBrewer_1.1-3    reticulate_1.47.0     parallelly_1.50.0    
    ##  [37] jquerylib_0.1.4       scattermore_1.2       Rcpp_1.1.2           
    ##  [40] knitr_1.52            future.apply_1.21.0   IRanges_2.46.0       
    ##  [43] sctransform_0.4.3     Matrix_1.7-6          igraph_2.3.4         
    ##  [46] tidyselect_1.2.1      abind_1.4-8           yaml_2.3.12          
    ##  [49] codetools_0.2-20      listenv_1.1.0         lattice_0.23-1       
    ##  [52] tibble_3.3.1          plyr_1.8.9            S7_1.0.0             
    ##  [55] evaluate_1.0.5        Rtsne_0.17            future_1.82.0        
    ##  [58] fastDummies_1.7.6     desc_1.4.3            pillar_1.11.1        
    ##  [61] MatrixGenerics_1.24.0 Seurat_6.0.0          KernSmooth_2.23-27   
    ##  [64] stats4_4.6.2          generics_0.1.4        RcppHNSW_0.7.0       
    ##  [67] sp_2.2-4              S4Vectors_0.50.3      ggplot2_4.0.4        
    ##  [70] scales_1.4.0          globals_0.20.0        glue_1.8.1           
    ##  [73] tools_4.6.2           data.table_1.18.6.1   RSpectra_0.16-2      
    ##  [76] SigBridgeRUtils_0.2.7 RANN_2.6.3            fs_2.1.0             
    ##  [79] dotCall64_1.2         cowplot_1.2.0         grid_4.6.2           
    ##  [82] patchwork_1.3.2       cli_3.6.6             textshaping_1.1.0    
    ##  [85] spam_2.11-4           S4Arrays_1.12.1       dplyr_1.2.1          
    ##  [88] uwot_0.2.5            gtable_0.3.6          pls_2.9-0            
    ##  [91] sass_0.5.0            digest_0.6.39         progressr_1.1.0      
    ##  [94] BiocGenerics_0.58.1   SparseArray_1.12.3    ggrepel_0.9.8        
    ##  [97] htmlwidgets_1.6.4     SeuratObject_5.4.0    farver_2.1.2         
    ## [100] htmltools_0.5.9       pkgdown_2.2.1         lifecycle_1.0.5      
    ## [103] httr_1.4.9            MASS_7.3-66

------------------------------------------------------------------------

### 6. References

> 1.  Sun D, Guan X, Moran AE, Wu LY, Qian DZ, Schedin P, et
>     al. Identifying phenotype-associated subpopulations by integrating
>     bulk and single-cell sequencing data. Nat Biotechnol. 2022
>     Apr;40(4):527-38.
>
> 2.  Xie A, Wang H, Zhao J, Wang Z, Xu J, Xu Y. scPAS: single-cell
>     phenotype-associated subpopulation identifier. Briefings in
>     Bioinformatics. 2024 Nov 22;26(1):bbae655.
>
> 3.  Zhang Q, Jin S, Zou X. scAB detects multiresolution cell states
>     with clinical significance by integrating single-cell genomics and
>     bulk sequencing data. Nucleic Acids Research. 2022 Nov
>     28;50(21):12112-30.
>
> 4.  He Y, Long R, Wang X. Inferring Phenotypes of Single Cells Based
>     on the Expression Profiles of Phenotype-Associated Marker Genes in
>     Bulks and Single Cells. Interdiscip Sci. Published online January
>     28, 2026.
>
> 5.  Johnson TS, Yu CY, Huang Z, Xu S, Wang T, Dong C, et
>     al. Diagnostic Evidence GAuge of Single cells (DEGAS): a flexible
>     deep transfer learning framework for prioritizing cells in
>     relation to disease. Genome Med. 2022 Feb 1;14(1):11.
>
> 6.  Li J, Zhang H, Mu B, Zuo H, Zhou K. Identifying
>     phenotype-associated subpopulations through LP_SGL. Briefings in
>     Bioinformatics. 2023 Nov 22;25(1):bbad424.
>
> 7.  Ruan X, Cheng Y, Ye Y, Wang Y, Chen X, Yang Y, et al. PIPET:
>     predicting relevant subpopulations in single-cell data using
>     phenotypic information from bulk data. Briefings in
>     Bioinformatics. 2024 May 23;25(4):bbae260.
>
> 8.  Jolasun, Y. et al. SIDISH integrates single-cell and bulk
>     transcriptomics to identify high-risk cells and guide precision
>     therapeutics through in silico perturbation. Nat Commun 16, 11271
>     (2025).
>
> 9.  Gan, D., Zhu, Y., Lu, X. & Li, J. SCIPAC: Quantitative estimation
>     of cell-phenotype associations. Genome Biol 25, 119 (2024).
>
> 10. Lin, Y. et al. TiRank prioritizes phenotypic niches in tumor
>     microenvironment for clinical biomarker discovery. Genome Med 18,
>     23 (2026).
