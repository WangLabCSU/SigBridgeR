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
