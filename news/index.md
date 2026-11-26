# Changelog

## SigBridgeR 4.0.0

### NEW FEATURES

- Added
  [`GetGeneLength()`](https://wanglabcsu.github.io/sigbridger/reference/GetGeneLength.md)
  and
  [`FPKMToTPM()`](https://wanglabcsu.github.io/sigbridger/reference/FPKMToTPM.md)
  for extracting gene lengths from GTF annotations and converting FPKM
  to TPM.
  [`GetGeneLength()`](https://wanglabcsu.github.io/sigbridger/reference/GetGeneLength.md)
  supports a `canonical_chr_only` option and preserves gene ID versions.
- Added
  [`DoPhenoMapR()`](https://wanglabcsu.github.io/sigbridger/reference/DoPhenoMapR.md)
  as a screening method for PhenoMapR, mapping bulk phenotypes to
  single-cell subpopulations, with
  [`ValidatePhenoMapRParams()`](https://wanglabcsu.github.io/sigbridger/reference/ValidatePhenoMapRParams.md)
  for input alignment and parameter checks. It is registered in
  `ScreenStrategy` and available through
  [`Screen()`](https://wanglabcsu.github.io/sigbridger/reference/Screen.md).
- Introduced a benchmarking framework:
  [`BenchAcc()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAcc.md),
  [`BenchGSEA()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEA.md),
  [`BenchGeneScore()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScore.md)
  and
  [`BenchMemory()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemory.md),
  along with their S7 configuration classes (`BenchConfig`,
  `BenchAccConfig`, `BenchGSEAConfig`, `BenchGeneScoreConfig`,
  `BenchMemoryConfig`, `MethodParamSet` and the `*ConfigSet` helpers).
  [`BenchAcc()`](https://wanglabcsu.github.io/sigbridger/reference/BenchAcc.md)
  now also reports `Precision` and `Recall`.
- Exported `AVAILABLE_SCREEN_METHODS` and `AVAILABLE_ANNOTATE_METHODS`;
  the `ScreenStrategy` registry replaces ad-hoc `names(ScreenStrategy)`
  lookups.
- [`SymbolConvert()`](https://wanglabcsu.github.io/sigbridger/reference/SymbolConvert.md)
  now accepts character vectors, offers `update_symbol_method`
  (`"scCustomize"`, `"Seurat"` or `"auto"`), supports extended genome
  builds, and replaces unmatched/duplicated symbols with the original
  IDs; `unknown_format` is deprecated and ignored.
- Added a `+` operator for merging `ScreenMethodResult` objects and S3
  methods for merging `Seurat` objects.
- [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  interactively prompts for `group.by` when multiple candidates are
  available.
- Added
  [`list_to_chr_formated()`](https://wanglabcsu.github.io/sigbridger/reference/list_to_chr_formated.md),
  which renders a list as an indented multi-line string;
  [`list_to_chr()`](https://wanglabcsu.github.io/sigbridger/reference/list_to_chr.md)
  and the [`print()`](https://rdrr.io/r/base/print.html) method for
  `SigBridgeRBase` objects are now backed by C++.
- [`SetThreads()`](https://wanglabcsu.github.io/sigbridger/reference/SetThreads.md)
  (renamed from
  [`setThreads()`](https://satijalab.org/seurat/reference/setThreads.html))
  gains a Seurat thread-config backend; the WGCNA backend is renamed
  `wgcna`.

### PERFORMANCE

- [`FindRobustElbow()`](https://wanglabcsu.github.io/sigbridger/reference/FindRobustElbow.md)
  offloads elbow detection to C++ using tatami streaming and caps
  `ndims` at the number of available PCs.
- [`FPKMToTPM()`](https://wanglabcsu.github.io/sigbridger/reference/FPKMToTPM.md)
  /
  [`GetGeneLength()`](https://wanglabcsu.github.io/sigbridger/reference/GetGeneLength.md)
  stream GTF files through a zlib buffered reader with single-pass C++
  validation.
- Replaced
  [`utils::modifyList()`](https://rdrr.io/r/utils/modifyList.html) with
  a custom C++
  [`modify_list()`](https://wanglabcsu.github.io/sigbridger/reference/modify_list.md)
  and added an `nthreads` option.
- [`BenchMemory()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemory.md)
  serialises its staged inputs with `compress_level = 1` and derives the
  `qs2` thread count once; it uses at least two cores and skips the
  machine-resource probe when `probe_resources = FALSE`.
- Matrix handling now uses `beachmat`.

### BREAKING CHANGES

- Minimum requirements raised: R \>= 4.2.0, S7 \>= 1.0.0, C++17 and zlib
  (\>= 1.2.4) (previously R \>= 4.1.0, S7 \>= 1.0.0, C++14).
- [`CacheSysCall()`](https://wanglabcsu.github.io/sigbridger/reference/CacheSysCall.md)
  load mode now requires an explicit `path`.
- [`SymbolConvert()`](https://wanglabcsu.github.io/sigbridger/reference/SymbolConvert.md)’s
  `unknown_format` argument is deprecated and ignored.
- `PIPET_Screen()` now defaults `distance` to `"cosine"`;
  [`DoDEGAS()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGAS.md)
  defaults `normality_test_method` to `"jarque-bera"`.
- [`setThreads()`](https://satijalab.org/seurat/reference/setThreads.html)
  is renamed
  [`SetThreads()`](https://wanglabcsu.github.io/sigbridger/reference/SetThreads.md);
  the WGCNA backend is renamed `wgcna`.
- The WGCNA backend now requires at least two threads: values below `2L`
  are raised to `2L` (with a warning when `verbose = TRUE`).

### BUG FIXES

- Fixed a Scissor crash (“a character vector argument expected”) caused
  by missing null-checks on cache paths
  ([\#22](https://github.com/WangLabCSU/SigBridgeR/issues/22)).
- Fixed bugs and integer-literal conversions in the vendored SCIPAC
  code.
- Added typed error messages across modules and switched to
  [`rlang::arg_match0()`](https://rlang.r-lib.org/reference/arg_match.html)
  /
  [`rlang::arg_match()`](https://rlang.r-lib.org/reference/arg_match.html)
  for stricter argument validation.
- Improved
  [`SymbolConvert()`](https://wanglabcsu.github.io/sigbridger/reference/SymbolConvert.md)
  handling of NA and duplicated symbols.
- [`MergeResult()`](https://wanglabcsu.github.io/sigbridger/reference/MergeResult.md)
  now handles duplicate columns and slot merging correctly.
- Improved error handling in
  [`ValidateDEGASv2Params()`](https://wanglabcsu.github.io/sigbridger/reference/ValidateDEGASv2Params.md)
  and
  [`ValidateDEGASParams()`](https://wanglabcsu.github.io/sigbridger/reference/ValidateDEGASParams.md).
- Fixed `ChooseNormalizationCheck()` layer detection and tightened
  weight validation.
- Corrected the `ssGSEA` dependency check and added `AUCell`, `GSVA`,
  `UCell` to Suggests.
- Fixed typos `rTIRank` -\> `rTiRank` and `qs2` -\> `qd_save`, and added
  a path-normalization helper.
- Corrected integer defaults for `z_score_cutoff` in
  [`DoPhenoMapR()`](https://wanglabcsu.github.io/sigbridger/reference/DoPhenoMapR.md)
  and for `n_cores` in
  [`BenchMemory()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemory.md).
- Removed a duplicate
  [`cli::cli_end()`](https://cli.r-lib.org/reference/cli_end.html) call
  in `S3_base`.
- Bumped the vendored Scissor fork to v2.1.0 with `qs2` cache support.
- Dropped the `randomcoloR` dependency and resolved R CMD check notes.

### MINOR IMPROVEMENTS

- Restructured the benchmark config hierarchy: shared properties moved
  to `BenchConfig`, with dedicated `BenchAccConfig` / `BenchGSEAConfig`
  classes.
- [`BenchMemoryConfigSet()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemoryConfigSet.md)
  gains a `probe_resources` flag to skip the machine-resource probe (and
  its [`gc()`](https://rdrr.io/r/base/gc.html)) before staging inputs.
- Added `S7_on_build()` / `S7_on_unload()` hooks, benchmark properties
  on `SigBridgeRBase`, and fully qualified S7 class references.
- Prefixed `rlang` calls with `rlang::` and replaced
  `cli::cli_fmt(cli::cli_text(...))` with
  [`cli::format_inline()`](https://cli.r-lib.org/reference/format_inline.html).
- Updated roxygen2 to 8.1.1 and added dynamic screen-method
  documentation generation.
- Added runtime `glue` dependency checks and raised the minimum `cli`
  version to 3.0.0.
- Added full third-party copyright attribution to `DESCRIPTION` and
  updated `CITATION`.
- Simplified named `...` validation in
  [`BenchGSEA()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGSEA.md)
  and
  [`BenchGeneScore()`](https://wanglabcsu.github.io/sigbridger/reference/BenchGeneScore.md)
  (`!all(nzchar(...))`).

### DOCUMENTATION

- Updated the README algorithm table and citation information.
- Refreshed roxygen documentation and benchmark docs covering C++
  backends and config APIs; reworked
  [`BenchMemory()`](https://wanglabcsu.github.io/sigbridger/reference/BenchMemory.md)
  docs for child-process measurement.
- Added benchmark vignettes for accuracy, GSEA, gene-score and memory
  benchmarking, and expanded `Full_Tutorial` and
  `Other_Function_Details`.
- Added
  [`DoPhenoMapR()`](https://wanglabcsu.github.io/sigbridger/reference/DoPhenoMapR.md)
  documentation and registered it in the algorithm tables.
- Added a spelling check to CI.
- Deprecated the `Extending_Seurat` vignette, which now carries a
  lifecycle badge and points to
  [`Register()`](https://wanglabcsu.github.io/sigbridger/reference/Register.md)
  and `SCIntegrate(pipeline = )`.

### INTERNAL CHANGES

- Added unit tests for benchmark functions and strategy registries, plus
  screen-method and
  [`AddMetaFeature()`](https://wanglabcsu.github.io/sigbridger/reference/AddMetaFeature.md)
  tests with fixtures.
- Added tests for
  [`DoPhenoMapR()`](https://wanglabcsu.github.io/sigbridger/reference/DoPhenoMapR.md)
  /
  [`ValidatePhenoMapRParams()`](https://wanglabcsu.github.io/sigbridger/reference/ValidatePhenoMapRParams.md),
  the C++
  [`list_to_chr()`](https://wanglabcsu.github.io/sigbridger/reference/list_to_chr.md)
  helpers, and the roxygen argument-documentation generator.
- CI now derives R-CMD-check dependencies from `DESCRIPTION`; updated
  the version-bumper workflow with cached
  [`cli::symbol`](https://cli.r-lib.org/reference/symbol.html).
- Added an AI skill.

## SigBridgeR 4.0.0-beta-4

### NEW FEATURES

- Added screening methods for DEGASv2 and TiRank.
- [`SymbolConvert()`](https://wanglabcsu.github.io/sigbridger/reference/SymbolConvert.md)
  gains an `update_symbol_method` argument to choose between
  `"scCustomize"` and `"Seurat"` for updating gene symbols (`"auto"`
  prefers `scCustomize` when installed).
- [`SymbolConvert()`](https://wanglabcsu.github.io/sigbridger/reference/SymbolConvert.md)
  now replaces unmatched IDs (NA gene symbols) with the original
  Ensembl/TCGA IDs; the `unknown_format` argument is deprecated and
  ignored.

### BUG FIXES

- Corrected DEGASv2 preprocessing and validation, including handling of
  zero-gene selections.

## SigBridgeR 4.0.0-beta-3-release

### NEW FEATURES

- **CellTypist integration** 🔬  
  CellTypist cell-type annotation is now powered by a native Python
  module.  
  New helper functions support:
  - downloading and listing available CellTypist models
  - selecting the assay to use for annotation  
    This makes running CellTypist from R seamless and more flexible.
- **[`summary()`](https://rdrr.io/r/base/summary.html) method for
  SigBridgeR objects**  
  Added a [`summary()`](https://rdrr.io/r/base/summary.html) S3 method
  that prints a concise overview of key metrics, screening results, and
  object metadata, improving data exploration.
- **Modular backends in
  [`setThreads()`](https://satijalab.org/seurat/reference/setThreads.html)
  and WGCNA support**  
  The
  [`setThreads()`](https://satijalab.org/seurat/reference/setThreads.html)
  function now uses a modular backend design.  
  It includes explicit support for WGCNA parallelisation, enabling
  consistent multi-threaded behaviour across different algorithms.
- **Overhauled class system and S3 methods**  
  The core object classes have been redesigned for long-term
  maintainability and extensibility. All public S3 methods have been
  consolidated and new ones introduced where appropriate. This work sets
  the stage for future v4.0.0 developments.

### PERFORMANCE

- **C++ acceleration for metric extraction and normalisation** ⚡  
  Metric extraction and normalisation, the computational bottlenecks in
  large-scale screenings, have been offloaded to C++. Users can expect
  dramatic speed improvements on datasets with thousands of features.

### BREAKING CHANGES

- **Class system rework**  
  The overhaul of internal object classes and S3 methods may break code
  that directly manipulates the internal structure of `SigBridgeR`
  objects.  
  Exported functions remain stable, but custom scripts accessing
  unexported S3 methods or object slots should be reviewed and adapted.

### BUG FIXES

- Fixed all R CMD check errors and warnings, including documentation
  titles, parameter descriptions, and example code.
- Resolved compatibility issues with updated dependencies after a
  dependency cleanup.
- Vignette setup chunks now use `eval = FALSE` to prevent runtime
  failures during package installation and build.
- Improved error handling and documentation in
  [`mLLMCelltypeAnnotate()`](https://wanglabcsu.github.io/sigbridger/reference/mLLMCelltypeAnnotate.md).

### MINOR IMPROVEMENTS

- Simplified
  [`AddMetaFeature()`](https://wanglabcsu.github.io/sigbridger/reference/AddMetaFeature.md)
  and refreshed its documentation.
- Refined strategy inspection and DEGAS logic for more accurate
  screening guidance.
- Updated `data.table` dependency and import statements to the latest
  compatible versions.
- Renamed internal utility `as.character.list` to `list_to_chr` to avoid
  potential S3 dispatch conflicts.
- Removed trailing commas from function arguments (code style).

### DOCUMENTATION

- Updated **README** and added a **CITATION** file so users can properly
  cite the package.
- Converted roxygen2 documentation from LaTeX `\itemize`/`\describe` to
  markdown list syntax, improving readability and pkgdown rendering.
- Modernised pkgdown configuration, vignettes, and internal keyword
  usage.

### INTERNAL CHANGES

- CI workflows now install additional system dependencies and R packages
  (`testthat`, `pak`, `remotes`, `tidydr`, `ggupset`) for more reliable
  R CMD check runs.
- Refactored threading backends to a modular design, simplifying future
  extensions.
- Cleaned up `DESCRIPTION` and `NAMESPACE` dependencies, resolving
  several compatibility glitches.
- Various internal reorganisations tied to the class system and S3
  method overhaul.

## SigBridgeR 4.0.0-beta-2

### DOCUMENTATION

- Updated README with improved citation information. 📖

### INTERNAL CHANGES

- Refactored project dependencies and validation logic for better
  maintainability.
- Performed general code cleanup and documentation updates.

## SigBridgeR 4.0.0-beta-1

### BREAKING CHANGES

- The computational core has been rewritten in C++ (via Rcpp), and the
  public API has been substantially restructured. Function and method
  signatures, exported operators, and class hierarchies have changed;
  review your scripts for compatibility.

### PERFORMANCE

- Core computations now run in compiled C++ code, yielding significant
  speed improvements.

### MINOR IMPROVEMENTS

- S7 class version information is now stored in a common base class,
  with improved validation.
- Error handling has been refined to provide clearer, more informative
  messages.

### DOCUMENTATION

- Vignettes and README are now distributed as plain Markdown files.
- Documentation has been updated to reflect the redesigned API.
- Chinese translations of documentation have been added.

### INTERNAL CHANGES

- Internal S7 class structure was refactored for better maintainability.
- Package dependencies were updated.
- The pkgdown site configuration and workflow were modernized.
- General code clean-up and routine maintenance.
