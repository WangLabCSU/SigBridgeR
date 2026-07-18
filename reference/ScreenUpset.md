# ScreenUpset - Visualize cell type intersections from screened Seurat object

This function creates an UpSet plot to visualize intersections between
different screening methods (e.g., scissor, scAB, scPAS, scPP) in a
Seurat object metadata. It calculates all possible combinations of
screening types and displays the number of cells positive for each
combination.

## Usage

``` r
ScreenUpset(
  screened_seurat,
  screen_type = NULL,
  order_by = c("freq", "degree"),
  show_plot = FALSE,
  n_intersections = 20,
  bar_color = "#4E79A7",
  combmatrix_point_color = "black",
  verbose = SigBridgeRUtils::getFuncOption("verbose"),
  ...
)
```

## Arguments

- screened_seurat:

  A Seurat object containing screening results in metadata. Must contain
  columns with screening types marked as "Positive".

- screen_type:

  Character vector of screening types to analyze. Default: NULL,
  indicating that a self-search pattern will be used.

- order_by:

  Character vector of columns to order the intersections by. Default:
  c("freq", "degree").

- show_plot:

  Whether to show the upset plot. Default: FALSE

- n_intersections:

  Number of intersections to display in the plot. Default: 20.

- bar_color:

  Color for the bars in the plot. Default: "#4E79A7".

- combmatrix_point_color:

  Color for points in the combination matrix. Default: "black".

- verbose:

  Logical, whether to print a message

- ...:

  Additional arguments passed to `ggupset::scale_x_upset()`,
  `ggupset::theme_combmatrix()`, and
  [`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html)
  for customizing the plot appearance. Argumetns are auto-filtered

## Value

A list containing two elements:

- `plot`: A ggplot object displaying the UpSet plot

- `stats`: A data frame with intersection statistics including:

  - `intersection`: Name of the intersection

  - `sets`: List of screening types in the intersection

  - `count`: Number of cells positive for all screening types in the
    intersection

## Details

The function performs the following steps:

1.  Validates input parameters and checks if specified screening types
    exist in metadata

2.  Generates all possible combinations of screening types (from 1 to
    the total number of types)

3.  Creates a logical matrix of positive cells for efficient computation

4.  Calculates cell counts for each intersection using vectorized
    operations

5.  Creates an UpSet plot visualization with customizable appearance

## See also

Other visualization_function:
[`ScreenFractionPlot()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenFractionPlot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with default parameters
result <- ScreenUpset(screened_seurat = my_seurat_obj)

# Customize screening types and appearance
result <- ScreenUpset(
  screened_seurat = my_seurat_obj,
  screen_type = c("scissor", "scAB"),
  n_intersections = 15,
  title = "Custom Title",
  bar_color = "darkred"
)
} # }

```
