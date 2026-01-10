# convert regex patterns to column names (internal)

An internal utility function that converts regular expression patterns
used for quality control in single-cell RNA-seq analysis into
standardized column names for Seurat object metadata. This function
handles both single patterns and combined patterns with logical OR
operators.

## Usage

``` r
Pattern2Colname(pat)
```

## Arguments

- pat:

  A character string containing a regular expression pattern or multiple
  patterns combined with `|` (OR operator).

## Value

A character string with the standardized column name derived from the
input pattern(s). The output follows these conventions: - Lowercase
letters only - Special characters removed or replaced with underscores -
Common patterns mapped to standardized abbreviations - Combined patterns
sorted alphabetically and joined with underscores

## See also

Other single_cell_preprocess:
[`FilterTumorCell()`](https://wanglabcsu.github.io/sigbridger/reference/FilterTumorCell.md),
[`FindRobustElbow()`](https://wanglabcsu.github.io/sigbridger/reference/FindRobustElbow.md),
[`QCPatternDetect()`](https://wanglabcsu.github.io/sigbridger/reference/QCPatternDetect.md),
[`SCPreProcess()`](https://wanglabcsu.github.io/sigbridger/reference/SCPreProcess.md),
[`compatible_with_3.0.2()`](https://wanglabcsu.github.io/sigbridger/reference/compatible_with_3.0.2.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Internal usage examples
Pattern2Colname("^MT-")                    # Returns "mt"
Pattern2Colname("^RP[LS]")                 # Returns "rp"
Pattern2Colname("^[rt]rna")                # Returns "rt_rna"
Pattern2Colname("^MT-|^RP[LS]")            # Returns "mt_rp"
Pattern2Colname("^HB[AB]?")                # Returns "hb_ab"
Pattern2Colname("Custom_Pattern[0-9]+")    # Returns "custom_pattern_0_9"
} # }
```
