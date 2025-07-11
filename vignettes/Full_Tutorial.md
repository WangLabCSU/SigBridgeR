# **Full Tutorial for SigBridgeR** {#full-tutorial-for-sigbridger}

## 0. Preface

### 0.1 Contents

- [**Full Tutorial for SigBridgeR** {#full-tutorial-for-sigbridger}](#full-tutorial-for-sigbridger-full-tutorial-for-sigbridger)
  - [0. Preface](#0-preface)
    - [0.1 Contents](#01-contents)
    - [0.1 Introduction to SigBridgeR](#01-introduction-to-sigbridger)
  - [1. Installation](#1-installation)
    - [1.1 Stable release from CRAN](#11-stable-release-from-cran)
    - [1.2 Development version from GitHub](#12-development-version-from-github)
    - [1.3 Check Dependencies](#13-check-dependencies)
  - [2. Loading and preprocessing data](#2-loading-and-preprocessing-data)
    - [2.1 Single-cell RNA-seq data](#21-single-cell-rna-seq-data)
      - [2.2.1 (Option A) Start from Seurat object](#221-option-a-start-from-seurat-object)
      - [2.2.2 (Option B) Start from raw matrix](#222-option-b-start-from-raw-matrix)
      - [2.2.3 (Option C) Start from AnnDataR6 object](#223-option-c-start-from-anndatar6-object)
      - [2.2.4 Preprocessing](#224-preprocessing)
    - [2.2 Bulk expression data](#22-bulk-expression-data)
    - [2.3 Mutational signature data](#23-mutational-signature-data)
    - [2.4 Matching Samples](#24-matching-samples)
  - [3. Runing SigBridgeR](#3-runing-sigbridger)
    - [3.1 (Option A) Scissor Screening](#31-option-a-scissor-screening)
    - [3.2 (Option B) scPAS Screening](#32-option-b-scpas-screening)
    - [3.3 (Option C) scAB Screening](#33-option-c-scab-screening)
    - [3.4 (Option D) scPP Screening](#34-option-d-scpp-screening)
  - [4. Visualization](#4-visualization)
    - [4.1 UMAP for screening results](#41-umap-for-screening-results)
    - [4.2 Stack bar plot for screening results](#42-stack-bar-plot-for-screening-results)
  - [5. Example](#5-example)
  - [6. Other function details](#6-other-function-details)
  - [7. References](#7-references)

### 0.1 Introduction to SigBridgeR

SigBridgeR (short for Mutational **Sig**nature **Bridge** in **R**) is an R package for screening tumour cell highly associated with mutational signatures from single-cell RNA-seq, bulk expression and mutational signature phenotype data at pan-cancer level. It is based on the R package [sunduanchen/Scissor](https://github.com/sunduanchen/Scissor), [Qinran-Zhang/scAB](https://github.com/Qinran-Zhang/scAB/), [WangX-Lab/ScPP](https://github.com/WangX-Lab/ScPP) and [aiminXie/scPAS](https://github.com/aiminXie/scPAS).

------------------------------------------------------------------------

## 1. Installation

You can install **SigBridgeR** using the following options:

### 1.1 Stable release from CRAN

```{r install_from_cran}
install.packages("SigBridgeR")
```

or

### 1.2 Development version from GitHub

```{r install_from_github}
if(!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("WangLabCSU/SigBridgeR")
```

### 1.3 Check Dependencies

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

------------------------------------------------------------------------

## 2. Loading and preprocessing data  

At first load `SigBridgeR` package:

```{r load_package}
library(SigBridgeR)
```

### 2.1 Single-cell RNA-seq data

#### 2.2.1 (Option A) Start from Seurat object

#### 2.2.2 (Option B) Start from raw matrix

#### 2.2.3 (Option C) Start from AnnDataR6 object

#### 2.2.4 Preprocessing 

### 2.2 Bulk expression data

### 2.3 Mutational signature data

### 2.4 Matching Samples

The original phenotype and bulk expression data may contain non-identical sample sets. This processing step performs sample matching to retain only the intersecting samples.

Some key details of `MatchSample`'s parameters:

- `ms_signature`: A data frame of mutational signatures after preprocessing
- `bulk_data`: A data frame of bulk expression data, 
- `col_id`: A character or numeric value specifying the column.

```{r matching_samples}
match_result<-MatchSample(
  ms_signature = tcga_ms_sbs,
  bulk_data = tcga_exp_count,
  col_id = col_id
)
```

After matching, the specified column in the mutation signature data (i.e., phenotypic data) will be binarized from continuous values - all positive numbers are recoded as 1.

The function returns a list containing:

- `phenotype`: A named vector of binary mutational signature data
- `matched_bulk_data`: A data frame of bulk expression data after intersecting samples
- `ms_select`: The column name specified, i.e. a particular mutation signature 

------------------------------------------------------------------------

## 3. Runing SigBridgeR

The function **`Screen`** provide 4 different options for screening mutational signatures, These four algorithms come from the repositories mentioned earlier, and you can choose one of them to screen your cells.:

### 3.1 (Option A) Scissor Screening

```{r scissor_screening}
scissor_result = Screen(
  matched_bulk = matched_bulk,
  sc_data = sc_dataset, # A Seurat object after preprocessing
  phenotype = matched_phenotype_data,
  label_type = "TP53", # The filtering labels are stored in the `\@misc` 
  phenotype_class = "binary", # 
  screen_method = c("Scissor"),
  path2save_scissor_inputs = "Tmp/Scissor_inputs.RData" # Intermediate data
)
```


### 3.2 (Option B) scPAS Screening

```{r scPAS_screening}

```


### 3.3 (Option C) scAB Screening

```{r scAB_screening}

```


### 3.4 (Option D) scPP Screening

```{r scPP_screening}

```


## 4. Visualization 

After screening, you can use these two functions for plotting the results of screening, **`FetchUMAP`** and **`ScreenFractionPlot`**.

### 4.1 UMAP for screening results

Some key details of `FetchUMAP`'s parameters:

-   `seurat_obj`: A Seurat object after screening.
-   `group_by`: Used to specify the column of the meta.data in `seurat_obj`. The plot results will be grouped by this parameter. Pass to `Seurat::DimPlot`'s `group.by` parameter.
-   `plot_color`: Custom color palette (named vector format):
    -   Required names: "Positive", "Negative", "Neutral"
    -   Default: c("Neutral"="#CECECE", "Positive"="#ff3333", "Negative"="#386c9b")
-   `order`: The order of the groups. Pass to `Seurat::DimPlot`'s `order` parameter.

Suppose you have performed `scissor` algorithm screening on your Seurat object and wish to examine the distribution across different celltypes and patient, you may reference and use the following code:

```{r umap_exmaple}
FetchUAMP(
  seurat_obj = your_seurat_obj,
  group_by = c("celltype","patient","scissor"),
  plot_color = NULL,
  plot_show = TRUE,
  order = c(2, 1),
)
```

This function generates three UMAP plots (one for each group specified in `group_by`), stored in a list. When `plot_show = TRUE`, you will see a composite plot displaying all groups together.

### 4.2 Stack bar plot for screening results

Some key details of `ScreenFractionPlot`'s parameters:

-   `seurat_obj`: A Seurat object after screening.
-   `group_by`: Used to specify the column of the meta.data in `seurat_obj`. The plot results will be grouped by this parameter.
-   `screen_type`: Screening algorithm used before. (case-sensitive, e.g., "scissor" for Scissor results)
-   `show_null`: Logical whether to show groups with zero cells (default: FALSE).
-   `plot_color` Custom color palette (named vector format):
    -   Required names: "Positive", "Negative", "Neutral"
    -   Default: c("Neutral"="#CECECE", "Positive"="#ff3333", "Negative"="#386c9b")

Suppose you have already performed the `scPAS` algorithm screening on your Seurat object, and you want to check the proportion of positive cells across different patients. You can refer to and use the following code

```{r stack_bar_plot_example}
plot <- ScreenFractionPlot(
   screened_seurat = your_seurat_obj,
   group_by = "patient", 
   screen_type = "scPAS",
   plot_title = "scPAS Screening Results"
)
```

The order of the groups is determined by the proportion of **Positive** cells within each group.

Use `?ScreenFractionPlot` in R terminal to see more details.

------------------------------------------------------------------------

## 5. Example

Here we use the example data to demonstrate how to use the functions in `SigBridgeR` to screen mutational signatures.


```{r example}
library(SigBridgeR)


```

------------------------------------------------------------------------

## 6. Other function details

-   `AddMisc()` : Add miscellaneous information to the Seurat object. Support for adding multiple attributes to the `SeuratObject@misc` slot simultaneously.

```{r add_misc_example}
# basic usage
seurat_obj <- AddMisc(seurat_obj, "QC_stats" = qc_df)

# Auto-incrementing example when `cover` set to FALSE
seurat_obj <- AddMisc(seurat_obj, markers = markers1)
seurat_obj <- AddMisc(seurat_obj, markers = markers2, cover=FALSE)

# Add multiple attributes to the `SeuratObject@misc` slot simultaneously
seurat_obj <- AddMisc(seurat_obj, markers1 = markers1, markers2 = markers2)
```

Use `?AddMisc` in R terminal to see more details.

## 7. References