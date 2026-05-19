# Package index

## Data Preprocessing

Functions for preprocessing inputs.

- [`BulkPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/BulkPreProcess.md)
  : Bulk RNA-seq Data Preprocessing and Quality Control Function
- [`SCPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcess.md)
  : Single-Cell RNA-seq Preprocessing Pipeline
- [`QCFilter()`](https://wanglabcsu.github.io/sigbridger/reference/QCFilter.md)
  : Filter Seurat object cells by QC metrics
- [`QCPatternDetect()`](https://wanglabcsu.github.io/sigbridger/reference/QCPatternDetect.md)
  : Calculate Percentage of Features Matching Patterns
- [`Pattern2Colname()`](https://wanglabcsu.github.io/sigbridger/reference/Pattern2Colname.md)
  : convert regex patterns to column names (internal)
- [`SymbolConvert()`](https://wanglabcsu.github.io/sigbridger/reference/SymbolConvert.md)
  : Convert Ensembles Version IDs & TCGA Version IDs to Genes in Bulk
  Expression Data
- [`AggregateDupRows()`](https://wanglabcsu.github.io/sigbridger/reference/aggregate-dups.md)
  [`AggregateDupCols()`](https://wanglabcsu.github.io/sigbridger/reference/aggregate-dups.md)
  [`AggregateDups()`](https://wanglabcsu.github.io/sigbridger/reference/aggregate-dups.md)
  : Aggregate Rows or Columns with Duplicate Names
- [`SCIntegrate()`](https://wanglabcsu.github.io/sigbridger/reference/SCIntegrate.md)
  [`SCIntegrate.data.frame()`](https://wanglabcsu.github.io/sigbridger/reference/SCIntegrate.md)
  [`SCIntegrate.Matrix()`](https://wanglabcsu.github.io/sigbridger/reference/SCIntegrate.md)
  [`SCIntegrate.Seurat()`](https://wanglabcsu.github.io/sigbridger/reference/SCIntegrate.md)
  : Integrate Single-Cell Datasets
- [`CheckNA()`](https://wanglabcsu.github.io/sigbridger/reference/CheckNA.md)
  : Check for Missing Values (NA) in Data
- [`PhenoPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/PhenoPreProcess.md)
  : Preprocess Phenotype Data
- [`PhenoMap()`](https://wanglabcsu.github.io/sigbridger/reference/PhenoMap.md)
  : Map Phenotype Values Using Conditional Rules

## Result Integration & Visualization

Merge and visualize screening results across methods or conditions.

- [`MergeResult()`](https://wanglabcsu.github.io/sigbridger/reference/MergeResult.md)
  : Merge Multiple Screening Analysis Results
- [`WeightedVote()`](https://wanglabcsu.github.io/sigbridger/reference/WeightedVote.md)
  : Weighted Voting Aggregation for Multi-Voter Classification
- [`ScreenFractionPlot()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenFractionPlot.md)
  : Visualization of Cell Screening Fractions
- [`ScreenUpset()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenUpset.md)
  : ScreenUpset - Visualize cell type intersections from screened Seurat
  object

## Screening Methods

Built-in screening algorithms.

- [`Screen()`](https://wanglabcsu.github.io/sigbridger/reference/Screen.md)
  : Single-Cell Data Screening
- [`DoLP_SGL()`](https://wanglabcsu.github.io/sigbridger/reference/DoLP_SGL.md)
  : Perform LP-SGL Screening Analysis
- [`DoPIPET()`](https://wanglabcsu.github.io/sigbridger/reference/DoPIPET.md)
  : Perform PIPET Screening Analysis
- [`DoScissor()`](https://wanglabcsu.github.io/sigbridger/reference/DoScissor.md)
  : Perform Scissor Screening Analysis
- [`DoscAB()`](https://wanglabcsu.github.io/sigbridger/reference/DoscAB.md)
  : Perform scAB Screening Analysis
- [`DoscPAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPAS.md)
  : Perform scPAS Screening Analysis
- [`DoscPP()`](https://wanglabcsu.github.io/sigbridger/reference/DoscPP.md)
  : Perform scPP screening analysis
- [`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md)
  : Run DEGAS Analysis for Single-Cell and Bulk RNA-seq Data Integration
- [`DoSIDISH()`](https://wanglabcsu.github.io/sigbridger/reference/DoSIDISH.md)
  : Perform SIDISH Screening Analysis
- [`DoSCIPAC()`](https://wanglabcsu.github.io/sigbridger/reference/DoSCIPAC.md)
  : Screen Single-Cell Data Using SCIPAC Algorithm
- [`DoTiRank()`](https://wanglabcsu.github.io/sigbridger/reference/DoTiRank.md)
  : Perform TiRank Screening Analysis

## Seurat Object Utilities

- [`AddMetaFeature()`](https://wanglabcsu.github.io/sigbridger/reference/AddMetaFeature.md)
  : Add Gene-Level Metadata to Seurat Object (Vectorized, ...-based)
- [`AddMisc()`](https://wanglabcsu.github.io/sigbridger/reference/AddMisc.md)
  : Safely Add Miscellaneous Data to Seurat Object
- [`FindRobustElbow()`](https://wanglabcsu.github.io/sigbridger/reference/FindRobustElbow.md)
  : Automatically determine optimal PCA dimensions using multiple robust
  methods
- [`ChooseNormalization()`](https://wanglabcsu.github.io/sigbridger/reference/ChooseNormalization.md)
  : Data-Driven Selection of Single-Cell Normalization Methods

## Annotate Cell Types in Single-cell Datasets

Built-in annotationa algorithms.

- [`SCAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/SCAnnotate.md)
  : Unified Interface for Single-Cell Annotation Methods
- [`CellTypistAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/CellTypistAnnotate.md)
  : Annotate Cell Types Using CellTypist (Python Backend)
- [`SingleRAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/SingleRAnnotate.md)
  : Annotate Single-Cell Data Using SingleR
- [`mLLMCelltypeAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/mLLMCelltypeAnnotate.md)
  : Annotate Cell Types Using Multi-LLM Consensus Approach

## Reference Data & External Resources

Load reference datasets.

- [`LoadRefData()`](https://wanglabcsu.github.io/sigbridger/reference/LoadRefData.md)
  : Download & Load Reference Data
- [`get_var_value()`](https://wanglabcsu.github.io/sigbridger/reference/get_var_value.md)
  : Trace and compute the value of a variable defined inside a function

## Python Environment Management

Manage Python environments for integrated R/Python workflows.

- [`SetupPyEnv()`](https://wanglabcsu.github.io/sigbridger/reference/SetupPyEnv.md)
  : Create or Use Python Environment with Required Packages
- [`ListPyEnv()`](https://wanglabcsu.github.io/sigbridger/reference/ListPyEnv.md)
  : List Available Python Environments

## Extension Tools

- [`InterceptStrategy()`](https://wanglabcsu.github.io/sigbridger/reference/InterceptStrategy.md)
  : Inspect Registered Strategy Environments
- [`Register()`](https://wanglabcsu.github.io/sigbridger/reference/Register.md)
  : Unified Registration Interface for Strategy Methods
- [`RegisterScreenMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterScreenMethod.md)
  : Register a Custom Screening Method for Phenotype-Driven Analysis
- [`RegisterAnnoMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterAnnoMethod.md)
  : Register an Annotation Method into the Strategy Registry
- [`RegisterSeuratMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterSeuratMethod.md)
  : Register a Seurat Processing Strategy
- [`ValidateScreenFunc()`](https://wanglabcsu.github.io/sigbridger/reference/ValidateScreenFunc.md)
  : Validate Custom Screening Function Compliance
- [`TemplateScreenFunc()`](https://wanglabcsu.github.io/sigbridger/reference/TemplateScreenFunc.md)
  : Generate a Template Screening Function with Optional Roxygen2
  Documentation

## Registry

Extended registry for SigBridgeR

- [`SCPreProcessStrategy`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcessStrategy.md)
  : Preprocessing Strategy Registry for Single-Cell Workflows
- [`SCAnnotateStrategy`](https://wanglabcsu.github.io/sigbridger/reference/SCAnnotateStrategy.md)
  : Registry of Cell Type Annotation Methods
- [`ScreenStrategy`](https://wanglabcsu.github.io/sigbridger/reference/ScreenStrategy.md)
  : Registry of Phenotype-Associated Cell Screening Methods

## Cache Management

Manage cache files for SigBridgeR

- [`WriteCacheMeta()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCacheMeta.md)
  : Write Cache Metadata
- [`CheckCache()`](https://wanglabcsu.github.io/sigbridger/reference/CheckCache.md)
  : Check Cache Configuration Consistency

## Package Configuration

Get and set global options for SigBridgeR behavior.

- [`getFuncOption()`](https://wanglabcsu.github.io/sigbridger/reference/getFuncOption.md)
  : Configuration Functions for SigBridgeR Package
- [`setFuncOption()`](https://wanglabcsu.github.io/sigbridger/reference/setFuncOption.md)
  : Configuration Functions for SigBridgeR Package
- [`setThreads()`](https://wanglabcsu.github.io/sigbridger/reference/setThreads.md)
  : Configure Parallel Execution Backends
