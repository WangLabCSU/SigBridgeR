# Validate and Prepare DoDEGASv2 Parameters

Internal helper that handles package installation checks, input
validation, and default-value resolution for
[`DoDEGASv2()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGASv2.md).

## Usage

``` r
ValidateDEGASv2Params(
  matched_bulk,
  phenotype,
  label_type = NULL,
  phenotype_class = c("binary", "survival", "continuous"),
  top_fraction_pos = 0.2,
  sclab = NULL,
  bulk_hvg = TRUE,
  bulk_de = TRUE,
  sc_de = TRUE,
  add_genes = NULL,
  n_hvg = 250L,
  n_bulk_de = 250L,
  n_sc_de = 200L,
  padj.thresh = 0.05,
  only.pos = FALSE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  n_st_classes = NULL,
  loss_type = c("cross_entropy", "log_neg", "rank_loss"),
  transfer_type = c("Wasserstein", "MMD"),
  model_save_dir = "DEGASv2_res",
  lambda1 = 1L,
  lambda2 = 3L,
  lambda3 = 3L,
  tot_seeds = 10L,
  tot_iters = 300L,
  extract_embs = FALSE,
  random_feat = FALSE,
  random_perc = 0.8,
  early_stopping = FALSE,
  ...
)
```

## Arguments

- matched_bulk, phenotype, label_type, phenotype_class:

  Forwarded from
  [`DoDEGASv2()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGASv2.md).

- top_fraction_pos, sclab, bulk_hvg, bulk_de, sc_de, add_genes:

  Forwarded from
  [`DoDEGASv2()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGASv2.md).

- n_hvg, n_bulk_de, n_sc_de, padj.thresh, only.pos, min.pct,
  logfc.threshold:

  Forwarded from
  [`DoDEGASv2()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGASv2.md).

- n_st_classes, loss_type, transfer_type, model_save_dir:

  Forwarded from
  [`DoDEGASv2()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGASv2.md).

- lambda1, lambda2, lambda3, tot_seeds, tot_iters:

  Forwarded from
  [`DoDEGASv2()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGASv2.md).

- extract_embs, random_feat, random_perc, early_stopping:

  Forwarded from
  [`DoDEGASv2()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGASv2.md).

- ...:

  Additional dots forwarded from
  [`DoDEGASv2()`](https://wanglabcsu.github.io/sigbridger/reference/DoDEGASv2.md).

## Value

A named list with elements: `phenotype_class`, `loss_type`,
`transfer_type`, `verbose`, `seed`.
