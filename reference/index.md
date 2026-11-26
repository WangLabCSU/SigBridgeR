# Package index

## Data Preprocessing

Functions for preprocessing inputs.

- [`SCPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcess.md)
  **\[deprecated\]** : Single-Cell RNA-seq Preprocessing Pipeline
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
  [`SCIntegrate.matrix()`](https://wanglabcsu.github.io/sigbridger/reference/SCIntegrate.md)
  [`SCIntegrate.Seurat()`](https://wanglabcsu.github.io/sigbridger/reference/SCIntegrate.md)
  : Integrate Single-Cell Datasets
- [`CheckNA()`](https://wanglabcsu.github.io/sigbridger/reference/CheckNA.md)
  : Check for Missing Values (NA) in Data
- [`PhenoPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/PhenoPreProcess.md)
  **\[deprecated\]** : Preprocess Phenotype Data
- [`PhenoMap()`](https://wanglabcsu.github.io/sigbridger/reference/PhenoMap.md)
  : Map Phenotype Values Using Conditional Rules
- [`CountsToTPM()`](https://wanglabcsu.github.io/sigbridger/reference/CountsToTPM.md)
  : Convert bulk RNA-seq raw counts to TPM
- [`FPKMToTPM()`](https://wanglabcsu.github.io/sigbridger/reference/FPKMToTPM.md)
  : Convert bulk RNA-seq FPKM to TPM
- [`GetGeneLength()`](https://wanglabcsu.github.io/sigbridger/reference/GetGeneLength.md)
  : Compute Gene Lengths from a GTF File or GTF Lines
- [`IsCountsMatrix()`](https://wanglabcsu.github.io/sigbridger/reference/IsCountsMatrix.md)
  : Check whether a matrix looks like RNA-seq raw counts
- [`IsSkewedDynamic()`](https://wanglabcsu.github.io/sigbridger/reference/IsSkewedDynamic.md)
  : Test Whether the Proportion of a Target Value Is Highly Skewed

## Result Integration & Visualization

Merge and visualize screening results across methods or conditions.

- [`MergeResult()`](https://wanglabcsu.github.io/sigbridger/reference/MergeResult.md)
  : Merge Multiple Screening Analysis Results

- [`` `+`( ``*`<SigBridgeR::ScreenMethodResult>`*`)`](https://wanglabcsu.github.io/sigbridger/reference/plus-methods.md)
  [`` `+`( ``*`<Seurat>`*`)`](https://wanglabcsu.github.io/sigbridger/reference/plus-methods.md)
  :

  Merge two objects with `+`

- [`WeightedVote()`](https://wanglabcsu.github.io/sigbridger/reference/WeightedVote.md)
  : Weighted Voting Aggregation for Multi-Voter Classification

- [`DrawBulkPCA()`](https://wanglabcsu.github.io/sigbridger/reference/DrawBulkPCA.md)
  : Draw PCA Plot for Bulk Expression Data

- [`ScreenFractionPlot()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenFractionPlot.md)
  : Visualization of Cell Screening Fractions

- [`ScreenUpset()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenUpset.md)
  : ScreenUpset - Visualize cell type intersections from screened Seurat
  object

- [`reexports`](https://wanglabcsu.github.io/sigbridger/reference/reexports.md)
  [`autoplot`](https://wanglabcsu.github.io/sigbridger/reference/reexports.md)
  : Objects exported from other packages

- [`autoplot(`*`<Seurat>`*`)`](https://wanglabcsu.github.io/sigbridger/reference/autoplot-SigBridgeR.md)
  [`autoplot(`*`<SigBridgeR::ScreenMethodResult>`*`)`](https://wanglabcsu.github.io/sigbridger/reference/autoplot-SigBridgeR.md)
  : autoplot method for Seurat UMAP

- [`palette_SigBridgeR()`](https://wanglabcsu.github.io/sigbridger/reference/palette_SigBridgeR.md)
  : A Palette

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
- [`DoDEGASv2()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGASv2.md)
  : Screen Function Template
- [`DoSIDISH()`](https://wanglabcsu.github.io/sigbridger/reference/DoSIDISH.md)
  : Perform SIDISH Screening Analysis
- [`DoSCIPAC()`](https://wanglabcsu.github.io/sigbridger/reference/DoSCIPAC.md)
  : Screen Single-Cell Data Using SCIPAC Algorithm
- [`DoTiRank()`](https://wanglabcsu.github.io/sigbridger/reference/DoTiRank.md)
  : Perform TiRank Screening Analysis
- [`DoPhenoMapR()`](https://wanglabcsu.github.io/sigbridger/reference/DoPhenoMapR.md)
  : PhenoMapR: Bulk-to-Single-Cell Phenotype Mapping

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
- [`DownloadCellTypistModels()`](https://wanglabcsu.github.io/sigbridger/reference/DownloadCellTypistModels.md)
  : Download CellTypist Models
- [`ListCellTypistModels()`](https://wanglabcsu.github.io/sigbridger/reference/ListCellTypistModels.md)
  : List Available CellTypist Models

## Reference Data & External Resources

Load reference datasets.

- [`LoadRefData()`](https://wanglabcsu.github.io/sigbridger/reference/LoadRefData.md)
  : Download & Load Reference Data

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
  **\[deprecated\]** : Register a Custom Screening Method for
  Phenotype-Driven Analysis
- [`RegisterAnnoMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterAnnoMethod.md)
  **\[deprecated\]** : Register an Annotation Method into the Strategy
  Registry
- [`RegisterSeuratMethod()`](https://wanglabcsu.github.io/sigbridger/reference/RegisterSeuratMethod.md)
  **\[deprecated\]** : Register a Seurat Processing Strategy
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
  [`AVAILABLE_ANNOTATE_METHODS`](https://wanglabcsu.github.io/sigbridger/reference/SCAnnotateStrategy.md)
  : Registry of Cell Type Annotation Methods
- [`ScreenStrategy`](https://wanglabcsu.github.io/sigbridger/reference/ScreenStrategy.md)
  [`AVAILABLE_SCREEN_METHODS`](https://wanglabcsu.github.io/sigbridger/reference/ScreenStrategy.md)
  : Registry of Phenotype-Associated Cell Screening Methods

## Benchmark Design

- [`BenchAcc()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAcc.md)
  : Benchmark Screening Methods: Classification Accuracy
- [`BenchAccConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAccConfig.md)
  : BenchAccConfig: Accuracy Benchmark Configuration
- [`BenchConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchConfig.md)
  : BenchConfig: Base Configuration for Benchmarks
- [`BenchGSEA()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEA.md)
  : Benchmark Screening Methods: Gene Set Enrichment Analysis
- [`BenchGSEAConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEAConfig.md)
  : BenchGSEAConfig: GSEA Benchmark Configuration
- [`BenchGSEAConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEAConfigSet.md)
  : Build GSEA Benchmark Controls
- [`BenchGeneScore()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScore.md)
  : Single-sample Gene Set Scoring for Single-Cell Data
- [`BenchGeneScoreConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfig.md)
  : BenchGeneScoreConfig: Gene Set Scoring Benchmark Configuration
- [`BenchGeneScoreConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScoreConfigSet.md)
  : Build Gene Set Scoring Benchmark Controls
- [`BenchMemory()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemory.md)
  : Benchmark Screening Methods: Wall Time and Memory
- [`BenchMemoryConfig()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfig.md)
  : BenchMemoryConfig: Time and Memory Benchmark Configuration
- [`BenchMemoryConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfigSet.md)
  : Resolve Benchmark Controls
- [`MethodParamSet()`](https://wanglabcsu.github.io/sigbridger/reference/MethodParamSet.md)
  : A List of Method Parameters

## Cache Management

Manage cache files for SigBridgeR

- [`WriteCacheMeta()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCacheMeta.md)
  : Write Cache Metadata
- [`CheckCache()`](https://wanglabcsu.github.io/sigbridger/reference/CheckCache.md)
  : Check Cache Configuration Consistency
- [`ChooseCache()`](https://wanglabcsu.github.io/sigbridger/reference/ChooseCache.md)
  : Choose a Cache Directory
- [`WriteCache()`](https://wanglabcsu.github.io/sigbridger/reference/WriteCache.md)
  : Write an R Object to Cache
- [`LoadCache()`](https://wanglabcsu.github.io/sigbridger/reference/LoadCache.md)
  : Load a Cached R Object
- [`CacheSetHere()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSetHere.md)
  : Resolve and Set Up Cache Directory
- [`CacheSysCall()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSysCall.md)
  : System-Level Cache Load and Save

## Package Configuration

Get and set global options for SigBridgeR behavior.

- [`getFuncOption()`](https://wanglabcsu.github.io/sigbridger/reference/getFuncOption.md)
  : Configuration Functions for SigBridgeR Package
- [`setFuncOption()`](https://wanglabcsu.github.io/sigbridger/reference/setFuncOption.md)
  : Configuration Functions for SigBridgeR Package
- [`SetThreads()`](https://wanglabcsu.github.io/sigbridger/reference/SetThreads.md)
  : Configure Parallel Execution Backends

## OOP Sytem Design

OOP system design for SigBridgeR

- [`SigBridgeRBase()`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeRBase.md)
  : SigBridgeRBase: Abstract Base Class for SigBridgeR
- [`Method()`](https://wanglabcsu.github.io/sigbridger/reference/Method.md)
  : Method: Abstract Base Class for Analysis Methods
- [`ScreenMethod()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethod.md)
  : ScreenMethod: Screening Method Class for Phenotype-Associated Cell
  Detection
- [`ScreenMethodCache()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodCache.md)
  : ScreenMethodCache: Screening Method Cache Container
- [`ScreenMethodConfig()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodConfig.md)
  : ScreenMethodConfig: Screening Method Configuration
- [`ScreenMethodResult()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenMethodResult.md)
  : ScreenMethodResult: Screening Method Result Container
- [`AnnotationMethod()`](https://wanglabcsu.github.io/sigbridger/reference/AnnotationMethod.md)
  : AnnotationMethod: Annotation Method Class for Cell Type
  Identification
- [`property_sigbridger_verison`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-properties.md)
  [`property_verison`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-properties.md)
  [`property_chr`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-properties.md)
  [`property_phenotype_class`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-properties.md)
  [`property_fn`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-properties.md)
  [`property_list`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-properties.md)
  [`property_mapper_fn`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-properties.md)
  [`property_py_env`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-properties.md)
  [`property_data_list`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-properties.md)
  [`property_seurat`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-properties.md)
  [`property_phenotype`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-properties.md)
  [`property_matrix`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-properties.md)
  : SigBridgeR Base Class Properties
- [`print(`*`<SigBridgeR::SigBridgeRBase>`*`)`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-S3.md)
  [`summary(`*`<SigBridgeR::SigBridgeRBase>`*`)`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-S3.md)
  [`length(`*`<SigBridgeR::SigBridgeRBase>`*`)`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-S3.md)
  [`names(`*`<SigBridgeR::SigBridgeRBase>`*`)`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-S3.md)
  [`` `names<-`( ``*`<SigBridgeR::SigBridgeRBase>`*`)`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-S3.md)
  [`format(`*`<SigBridgeR::SigBridgeRBase>`*`)`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-S3.md)
  [`as.list(`*`<SigBridgeR::SigBridgeRBase>`*`)`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-S3.md)
  [`as.character(`*`<SigBridgeR::SigBridgeRBase>`*`)`](https://wanglabcsu.github.io/sigbridger/reference/SigBridgeR-S3.md)
  : Base R Methods for SigBridgeR Objects
- [`` `$`( ``*`<SigBridgeR::SigBridgeRBase>`*`)`](https://wanglabcsu.github.io/sigbridger/reference/Property-visitor.md)
  [`` `[[`( ``*`<SigBridgeR::SigBridgeRBase>`*`)`](https://wanglabcsu.github.io/sigbridger/reference/Property-visitor.md)
  : Retrieve a property from SigBridgeR object
- [`list_to_chr()`](https://wanglabcsu.github.io/sigbridger/reference/list_to_chr.md)
  : Convert a List to a Compact Character Representation
- [`list_to_chr_formated()`](https://wanglabcsu.github.io/sigbridger/reference/list_to_chr_formated.md)
  : Convert a List to an Indented Multi-line Character Representation
