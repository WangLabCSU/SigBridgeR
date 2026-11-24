# SigBridgeR 4.0.0

## NEW FEATURES

* Added `GetGeneLength()` and `FPKMToTPM()` for extracting gene lengths from GTF annotations and converting FPKM to TPM. `GetGeneLength()` supports a `canonical_chr_only` option and preserves gene ID versions.
* Introduced a benchmarking framework: `BenchAcc()`, `BenchGSEA()`, `BenchGeneScore()` and `BenchMemory()`, along with their S7 configuration classes (`BenchConfig`, `BenchAccConfig`, `BenchGSEAConfig`, `BenchGeneScoreConfig`, `BenchMemoryConfig`, `MethodParamSet` and the `*ConfigSet` helpers). `BenchAcc()` now also reports `Precision` and `Recall`.
* Exported `AVAILABLE_SCREEN_METHODS` and `AVAILABLE_ANNOTATE_METHODS`; the `ScreenStrategy` registry replaces ad-hoc `names(ScreenStrategy)` lookups.
* `SymbolConvert()` now accepts character vectors, offers `update_symbol_method` (`"scCustomize"`, `"Seurat"` or `"auto"`), supports extended genome builds, and replaces unmatched/duplicated symbols with the original IDs; `unknown_format` is deprecated and ignored.
* Added a `+` operator for merging `ScreenMethodResult` objects and S3 methods for merging `Seurat` objects.
* `autoplot()` interactively prompts for `group.by` when multiple candidates are available.
* `setThreads()` gains a Seurat thread-config backend; the WGCNA backend is renamed `wgcna`.

## PERFORMANCE

* `FindRobustElbow()` offloads elbow detection to C++ using tatami streaming and caps `ndims` at the number of available PCs.
* `FPKMToTPM()` / `GetGeneLength()` stream GTF files through a zlib buffered reader with single-pass C++ validation.
* Replaced `utils::modifyList()` with a custom C++ `modify_list()` and added an `nthreads` option.
* Matrix handling now uses `beachmat`.

## BREAKING CHANGES

* Minimum requirements raised: R >= 4.2.0, S7 >= 1.0.0, C++17 and zlib (>= 1.2.4) (previously R >= 4.1.0, S7 >= 1.0.0, C++14).
* `CacheSysCall()` load mode now requires an explicit `path`.
* `SymbolConvert()`'s `unknown_format` argument is deprecated and ignored.
* `PIPET_Screen()` now defaults `distance` to `"cosine"`; `DoDEGAS()` defaults `normality_test_method` to `"jarque-bera"`.
* The WGCNA backend in `setThreads()` is renamed `wgcna`.

## BUG FIXES

* Fixed a Scissor crash ("a character vector argument expected") caused by missing null-checks on cache paths (#22).
* Fixed bugs and integer-literal conversions in the vendored SCIPAC code.
* Added typed error messages across modules and switched to `rlang::arg_match0()` / `rlang::arg_match()` for stricter argument validation.
* Improved `SymbolConvert()` handling of NA and duplicated symbols.
* `MergeResult()` now handles duplicate columns and slot merging correctly.
* Improved error handling in `ValidateDEGASv2Params()` and `ValidateDEGASParams()`.
* Fixed `ChooseNormalizationCheck()` layer detection and tightened weight validation.
* Corrected the `ssGSEA` dependency check and added `AUCell`, `GSVA`, `UCell` to Suggests.
* Fixed typos `rTIRank` -> `rTiRank` and `qs2` -> `qd_save`, and added a path-normalization helper.
* Removed a duplicate `cli::cli_end()` call in `S3_base`.
* Bumped the vendored Scissor fork to v2.1.0 with `qs2` cache support.
* Dropped the `randomcoloR` dependency and resolved R CMD check notes.

## MINOR IMPROVEMENTS

* Restructured the benchmark config hierarchy: shared properties moved to `BenchConfig`, with dedicated `BenchAccConfig` / `BenchGSEAConfig` classes.
* Added `S7_on_build()` / `S7_on_unload()` hooks, benchmark properties on `SigBridgeRBase`, and fully qualified S7 class references.
* Prefixed `rlang` calls with `rlang::` and replaced `cli::cli_fmt(cli::cli_text(...))` with `cli::format_inline()`.
* Updated roxygen2 to 8.1.1 and added dynamic screen-method documentation generation.
* Added runtime `glue` dependency checks and raised the minimum `cli` version to 3.0.0.
* Added full third-party copyright attribution to `DESCRIPTION` and updated `CITATION`.

## DOCUMENTATION

* Updated the README algorithm table and citation information.
* Refreshed roxygen documentation and benchmark docs covering C++ backends and config APIs; reworked `BenchMemory()` docs for child-process measurement.
* Added a spelling check to CI.
* Deprecated the `Extending_Seurat` vignette, which now carries a lifecycle badge and points to `Register()` and `SCIntegrate(pipeline = )`.

## INTERNAL CHANGES

* Added unit tests for benchmark functions and strategy registries, plus screen-method and `AddMetaFeature()` tests with fixtures.
* CI now derives R-CMD-check dependencies from `DESCRIPTION`; updated the version-bumper workflow with cached `cli::symbol`.
* Updated `.gitignore` / `.Rbuildignore`, excluded `dev/` from the package build, and removed hardcoded R binary paths from VS Code settings.
* Added an AI skill.

# SigBridgeR 4.0.0-beta-4

## NEW FEATURES

* Added screening methods for DEGASv2 and TiRank.
* `SymbolConvert()` gains an `update_symbol_method` argument to choose between `"scCustomize"` and `"Seurat"` for updating gene symbols (`"auto"` prefers `scCustomize` when installed).
* `SymbolConvert()` now replaces unmatched IDs (NA gene symbols) with the original Ensembl/TCGA IDs; the `unknown_format` argument is deprecated and ignored.

## BUG FIXES

* Corrected DEGASv2 preprocessing and validation, including handling of zero-gene selections.

# SigBridgeR 4.0.0-beta-3-release

## NEW FEATURES

*   **CellTypist integration** 🔬  
    CellTypist cell-type annotation is now powered by a native Python module.  
    New helper functions support:
    - downloading and listing available CellTypist models
    - selecting the assay to use for annotation  
    This makes running CellTypist from R seamless and more flexible.
*   **`summary()` method for SigBridgeR objects**  
    Added a `summary()` S3 method that prints a concise overview of key metrics, 
    screening results, and object metadata, improving data exploration.
*   **Modular backends in `setThreads()` and WGCNA support**  
    The `setThreads()` function now uses a modular backend design.  
    It includes explicit support for WGCNA parallelisation, enabling consistent 
    multi-threaded behaviour across different algorithms.
*   **Overhauled class system and S3 methods**  
    The core object classes have been redesigned for long-term maintainability 
    and extensibility. All public S3 methods have been consolidated and new ones 
    introduced where appropriate. This work sets the stage for future v4.0.0 
    developments.

## PERFORMANCE

*   **C++ acceleration for metric extraction and normalisation** ⚡  
    Metric extraction and normalisation, the computational bottlenecks in large-scale 
    screenings, have been offloaded to C++. Users can expect dramatic speed 
    improvements on datasets with thousands of features.

## BREAKING CHANGES

*   **Class system rework**  
    The overhaul of internal object classes and S3 methods may break code that 
    directly manipulates the internal structure of `SigBridgeR` objects.  
    Exported functions remain stable, but custom scripts accessing unexported 
    S3 methods or object slots should be reviewed and adapted.

## BUG FIXES

*   Fixed all R CMD check errors and warnings, including documentation titles, 
    parameter descriptions, and example code.
*   Resolved compatibility issues with updated dependencies after a dependency 
    cleanup.
*   Vignette setup chunks now use `eval = FALSE` to prevent runtime failures 
    during package installation and build.
*   Improved error handling and documentation in `mLLMCelltypeAnnotate()`.

## MINOR IMPROVEMENTS

*   Simplified `AddMetaFeature()` and refreshed its documentation.
*   Refined strategy inspection and DEGAS logic for more accurate screening 
    guidance.
*   Updated `data.table` dependency and import statements to the latest 
    compatible versions.
*   Renamed internal utility `as.character.list` to `list_to_chr` to avoid 
    potential S3 dispatch conflicts.
*   Removed trailing commas from function arguments (code style).

## DOCUMENTATION

*   Updated **README** and added a **CITATION** file so users can properly cite 
    the package.
*   Converted roxygen2 documentation from LaTeX `\itemize`/`\describe` to 
    markdown list syntax, improving readability and pkgdown rendering.
*   Modernised pkgdown configuration, vignettes, and internal keyword usage.

## INTERNAL CHANGES

*   CI workflows now install additional system dependencies and R packages 
    (`testthat`, `pak`, `remotes`, `tidydr`, `ggupset`) for more reliable 
    R CMD check runs.
*   Refactored threading backends to a modular design, simplifying future 
    extensions.
*   Cleaned up `DESCRIPTION` and `NAMESPACE` dependencies, resolving several 
    compatibility glitches.
*   Various internal reorganisations tied to the class system and S3 method 
    overhaul.

# SigBridgeR 4.0.0-beta-2

## DOCUMENTATION

* Updated README with improved citation information. 📖

## INTERNAL CHANGES

* Refactored project dependencies and validation logic for better maintainability.
* Performed general code cleanup and documentation updates.

# SigBridgeR 4.0.0-beta-1

## BREAKING CHANGES
* The computational core has been rewritten in C++ (via Rcpp), and the public API has been substantially restructured. Function and method signatures, exported operators, and class hierarchies have changed; review your scripts for compatibility.

## PERFORMANCE
* Core computations now run in compiled C++ code, yielding significant speed improvements.

## MINOR IMPROVEMENTS
* S7 class version information is now stored in a common base class, with improved validation.
* Error handling has been refined to provide clearer, more informative messages.

## DOCUMENTATION
* Vignettes and README are now distributed as plain Markdown files.
* Documentation has been updated to reflect the redesigned API.
* Chinese translations of documentation have been added.

## INTERNAL CHANGES
* Internal S7 class structure was refactored for better maintainability.
* Package dependencies were updated.
* The pkgdown site configuration and workflow were modernized.
* General code clean-up and routine maintenance.
