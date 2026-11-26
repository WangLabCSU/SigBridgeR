# Visualize the candidates proposed by robust elbow detection

Internal helper used by
[`FindRobustElbow()`](https://wanglabcsu.github.io/sigbridger/reference/FindRobustElbow.md)
to draw the variance curve of the PCA dimensions, the per-method
candidate cut-offs, and their labels.

## Usage

``` r
plot_robust_elbow(res)
```

## Arguments

- res:

  A list as returned by
  [`find_robust_elbow_cpp()`](https://wanglabcsu.github.io/sigbridger/reference/find_robust_elbow_cpp.md),
  containing at least `pct_variance`, `cumulative_variance`,
  `cumulative_90`, `cumulative_80`, `above_mean`, `above_2sd`,
  `second_deriv` and `distance`.

## Value

A `ggplot2` object.
