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
