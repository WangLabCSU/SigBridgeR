# Visualization of Cell Screening Fractions

Generates stacked bar plots showing the proportion of cells classified
as Positive/Negative/Neutral by single-cell screening algorithms
(Scissor, scPAS, scPP, or scAB) across different sample groups. Supports
both single and multiple screen types visualization.

## Usage

``` r
ScreenFractionPlot(
  screened_seurat,
  group_by = "Source", # x-axis label
  screen_type = NULL, # default: search all available screens
  show_null = FALSE,
  plot_color = NULL, # stack bar color of each status
  show_plot = FALSE,
  plot_title = "Screen Fraction",
  stack_width = 0.85,
  x_text_angle = 45,
  axis_linewidth = 0.8,
  x_lab = NULL,
  y_lab = "Fraction of Status",
  ncol = 2, # number of columns for facet wrap
  nrow = NULL, # number of rows for facet wrap
  verbose = SigBridgeRUtils::getFuncOption("verbose"),
  ... # ggplot2 theme arguments
)
```

## Arguments

- screened_seurat:

  A Seurat object containing screening results in metadata. Must contain
  columns corresponding to `screen_type`.

- group_by:

  Metadata column name for grouping samples (default: "Source").

- screen_type:

  Screening algorithm(s) used. Can be a single value or vector. Must
  match metadata column(s) (case-sensitive, e.g., "scissor" for Scissor
  results).

- show_null:

  Logical whether to show groups with zero cells (default: FALSE).

- plot_color:

  Custom color palette (named vector format): - Required names:
  "Positive", "Negative", "Neutral" - Default: c("Neutral"="#CECECE",
  "Other"="#CECECE", Positive"="#ff3333", "Negative"="#386c9b")

- show_plot:

  Logical to immediately display plot (default: TRUE).

- plot_title:

  Plot title (default: "Screen Fraction"). When multiple screen types,
  can be a vector of titles or single title (will append screen type
  prefix).

- stack_width:

  Bar width (default: 0.85).

- x_text_angle:

  X-axis label angle (default: 45).

- axis_linewidth:

  Axis line thickness (default: 0.8).

- x_lab:

  X-axis label (default: NULL).

- y_lab:

  Y-axis label (default: "Fraction of Status").

- ncol:

  Number of columns for facet wrap when multiple screen types (default:
  2).

- nrow:

  Number of rows for facet wrap when multiple screen types (default:
  NULL).

- verbose:

  Logical, whether to print a message

- ...:

  Other arguments passed to ggplot2::theme()

## Value

A list containing:

- stats: A data frame (single screen) or list of data frames (multiple
  screens) with screening statistics including:

  - Grouping variable counts

  - Raw cell counts

  - Percentage fractions

- plot: A ggplot2 object (single screen) or list of ggplot2 objects
  (multiple screens)

- combined_plot: A combined plot using patchwork (only for multiple
  screens)

## Visualization Details

- Bars are ordered by descending Positive fraction

- Y-axis shows percentage (0-100%)

- Zero-fraction groups are automatically hidden unless `show_null=TRUE`

- For multiple screen types, plots can be combined using patchwork

## See also

Other visualization_function:
[`ScreenUpset()`](https://wanglabcsu.github.io/sigbridger/reference/ScreenUpset.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Single screen type usage
res <- ScreenFractionPlot(
  screened_seurat = scissor_result,
  group_by = "PatientID",
  screen_type = "scissor"
)

# Multiple screen types at once
multi_res <- ScreenFractionPlot(
  screened_seurat = multi_screened_result,
  group_by = "TissueType",
  screen_type = c("scissor", "scPAS", "scPP", "scAB", "DEGAS"),
  ncol = 2
)
} # }
```
