> It is recommended to read the [Quick
> Start](https://wanglabcsu.github.io/SigBridgeR/articles/Quick_Start.html)
> or [Full
> Tutorial](https://wanglabcsu.github.io/SigBridgeR/articles/Full_Tutorial.html)
> first, as this will help understand the workflow of SigBridgeR.

## Introduction

With so many parameters in the screening methods, it’s often confusing
how to select the best ones for your data. This document provides some
tips for parameter tuning.

## Find Optimal Parameters for `Scissor`

To enable scissor to search an optimal `alpha`, just change
`alpha = <A value>` to `alpha = NULL`. In this case, alpha will iterate
from 0.005 to 0.9 until the selected cells fraction is over the `cutoff`
value. That is to say, Scissor will designate exactly the top (and
bottom) percentage of cells —- corresponding to the specified cutoff
value -— as positive cells.

    # * From
    # scissor_res <- Screen(
    #   matched_bulk = matched_bulk,
    #   sc_data = sc_dataset,
    #   phenotype = matched_phenotype_data,
    #   phenotype_class = "survival",
    #   screen_method = c("Scissor"),
    #   alpha = 0.05
    # )

    # * To
    scissor_res <- Screen(
      matched_bulk = matched_bulk,
      sc_data = sc_dataset,
      phenotype = matched_phenotype_data,
      phenotype_class = "survival",
      screen_method = c("Scissor"),
      alpha = NULL # <<< HERE
    )

Output goes like this when using example data

    # ℹ [2025/11/07 08:39:59] Scissor start...
    # ℹ [2025/11/07 08:39:59] Start from raw data...
    # ℹ Using "RNA_snn" graph for network.
    # ℹ [2025/11/07 08:39:59] Normalizing quantiles of data...
    # ℹ [2025/11/07 08:40:00] Subsetting data...
    # ℹ [2025/11/07 08:40:01] Calculating correlation...
    # --------------------------------------------------------------------------------------------
    # Five-number summary of correlations:
    # 0.105058 0.17552 0.192775 0.211366 0.325537
    # --------------------------------------------------------------------------------------------
    # ℹ [2025/11/07 08:40:05] Perform cox regression on the given clinical outcomes...
    # ✔ [2025/11/07 08:40:06] Statistics data saved to Scissor_inputs.RData.
    # ℹ [2025/11/07 08:40:07] Screening...

    # ── At alpha = 0.005 ──

    # Scissor identified 550 Scissor+ cells and 424 Scissor- cells.
    # The percentage of selected cell is: 89.113%

    # ── At alpha = 0.01 ──

    # Scissor identified 512 Scissor+ cells and 374 Scissor- cells.
    # The percentage of selected cell is: 81.061%

    # ── At alpha = 0.05 ──

    # Scissor identified 265 Scissor+ cells and 73 Scissor- cells.
    # The percentage of selected cell is: 30.924%

    # ── At alpha = 0.1 ──

    # Scissor identified 182 Scissor+ cells and 2 Scissor- cells.
    # The percentage of selected cell is: 16.834%
    # ℹ [2025/11/07 08:40:39] Scissor Ended.

## Find Optimal Parameters for `scAB`

scAB has two key parameters, `alpha` and `alpha_2`: the former
determines the overall selection proportion, while the latter determines
the selection proportion within each scAB subset. Simply change the
values of these two parameters from single values to numeric sequences,
and scAB will automatically iterate through the combinations to find the
optimal pair of `alpha` and `alpha_2`.

    # * From
    # scAB_res <- Screen(
    #   matched_bulk = matched_bulk,
    #   sc_data = sc_dataset,
    #   phenotype = matched_phenotype_data,
    #   phenotype_class = "survival",
    #   screen_method = c("scAB"),
    #   alpha = 0.05,
    #   alpha_2 = 0.05
    # )

    # * To
    seq1 <- c(0.005, 0.01, 0.05)
    seq2 <- c(5e-7, 1e-6, 5e-6, 1e-5)

    scAB_res <- Screen(
      matched_bulk = matched_bulk,
      sc_data = sc_dataset,
      phenotype = matched_phenotype_data,
      phenotype_class = "survival",
      screen_method = c("scAB"),
      alpha = seq1, # <<< HERE
      alpha_2 = seq2 # <<< HERE
    )

    length(seq1) * length(seq2)
    # [1] 12

As you can see, with so many (alpha, alpha\_2) pairs – each requiring
evaluation -— the runtime can become prohibitively long. In such cases,
parallel computing can be employed to accelerate the process.

    # * install `furrr` first, as it is required for parallel computing
    # * We recommend using `future.mirai`, as it is more stable and faster than `future`
    rlang::check_installed(c("furrr", "future.mirai"))

    future::plan(future.mirai::mirai_multisession, workers = 2L)

    scAB_res <- Screen(
      matched_bulk = matched_bulk,
      sc_data = sc_dataset,
      phenotype = matched_phenotype_data,
      phenotype_class = "survival",
      screen_method = c("scAB"),
      alpha = seq1,
      alpha_2 = seq2
      # parallel = TRUE # Auto-enabled
    )

Output goes like this when using example data

    # ℹ [2025/11/07 08:43:48] Start scAB screening.
    # ℹ  Using "RNA_snn" graph for network.
    # ℹ [2025/11/07 08:43:53] Selecting K...
    # • loss of 2: 0.003226
    # • loss of 3: 0.003273
    # ℹ [2025/11/07 08:44:04] Run NMF with phenotype and cell-cell similarity regularization at K = 2
    # ℹ [2025/11/07 08:44:04] Selecting optimal `alpha` and `alpha_2`, this would take a while
    # ℹ [2025/11/07 08:44:07] Using parallel processing
    #  Progress: ────────────────────────────────────────────────────────────────────────────────────────── 100%
    # ✔ [2025/11/07 08:46:08] Best `alpha` and `alpha_2` are 0.05 and 1e-05
    # ℹ [2025/11/07 08:46:14] Screening cells...
    # ℹ [2025/11/07 08:46:14] scAB screening done.

After identifying the optimal (alpha, alpha\_2) pair, scAB will use this
parameter pair to perform the final cell screening.

## Find Optimal Parameters for `scPAS`

The iteration of alpha in scPAS is similar to that in Scissor, and the
method to enable iteration is the same.

    # * From
    # scPAS_res <- Screen(
    #   matched_bulk = matched_bulk,
    #   sc_data = sc_dataset,
    #   phenotype = matched_phenotype_data,
    #   phenotype_class = "survival",
    #   screen_method = c("scPAS"),
    #   alpha = 0.05
    # )

    # * To
    scPAS_res <- Screen(
      matched_bulk = matched_bulk,
      sc_data = sc_dataset,
      phenotype = matched_phenotype_data,
      phenotype_class = "survival",
      screen_method = c("scPAS"),
      alpha = NULL # <<< HERE
    )

Output goes like this when using example data

    # ℹ [2025/10/20 16:43:31] Start scPAS screening.
    # ℹ [2025/10/20 16:43:32] Quantile normalization of bulk data.
    # ℹ [2025/10/20 16:43:32] Extracting single-cell expression profiles...
    # ℹ [2025/10/20 16:43:32] Constructing a gene-gene similarity by single cell data...
    # Building SNN based on a provided distance matrix
    # Computing SNN
    # ℹ [2025/10/20 16:43:33] Optimizing the network-regularized sparse regression model...
    # ℹ [2025/10/20 16:43:33] Perform cox regression on the given phenotypes...
    #
    # ── At alpha = 0.001 ──
    #
    # lambda = 7.61825638357188
    # scPAS identified 315 risk+ features and 366 risk- features.
    # The percentage of selected feature is: 77.918%
    #
    # ── At alpha = 0.005 ──
    #
    # lambda = 3.68743152046651
    # scPAS identified 193 risk+ features and 231 risk- features.
    # The percentage of selected feature is: 48.513%
    #
    # ── At alpha = 0.01 ──
    #
    # lambda = 2.55333682907113
    # scPAS identified 137 risk+ features and 158 risk- features.
    # The percentage of selected feature is: 33.753%
    #
    # ── At alpha = 0.05 ──
    #
    # lambda = 0.776168989003421
    # scPAS identified 59 risk+ features and 61 risk- features.
    # The percentage of selected feature is: 13.73%
    # ℹ [2025/10/20 16:44:55] Calculating quantified risk scores...
    # ℹ [2025/10/20 16:44:55] Qualitative identification by permutation test program with 2000 times random perturbations...
    # ✔ [2025/10/20 16:44:57] scPAS screening done.

## Find Optimal Parameters for `scPP`

scPP is based on marker genes, AUC, and enrichment scores, which
requires tuning the `probs` parameter, as it represents the quantile
threshold. The input `probs` must be between 0 and 0.5.

    # * From
    # scPP_res <- Screen(
    #   matched_bulk = matched_bulk,
    #   sc_data = sc_dataset,
    #   phenotype = matched_phenotype_data,
    #   phenotype_class = "survival",
    #   screen_method = c("scPP"),
    #   probs = 0.2
    # )

    # * To
    # * Note that since a future strategy was previously enabled, parallel computation will also be used here
    scPP_res <- Screen(
      matched_bulk = matched_bulk,
      sc_data = sc_dataset,
      phenotype = matched_phenotype_data,
      phenotype_class = "survival",
      screen_method = c("scPP"),
      probs = NULL # <<< HERE
    )

Output goes like this when using example data

    # ℹ [2025/11/07 08:53:20] Start scPP screening.
    # ℹ [2025/11/07 08:53:20] Finding overall markers...
    # Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights,  :
    #   Loglik converged before variable  1 ; coefficient may be infinite.
    # ℹ [2025/11/07 08:53:43] Screening...
    # ℹ [2025/11/07 08:53:43] Running optimization mode: testing 6 thresholds
    # ℹ [2025/11/07 08:53:43] Computing AUC scores...
    # Genes in the gene sets NOT available in the dataset:
    #   gene_pos:   13 (6% of 230)
    #   gene_neg:   54 (12% of 446)
    # ℹ [2025/11/07 08:53:44] Testing thresholds and computing NES differences...
    # ℹ [2025/11/07 08:53:46] Using parallel processing
    #  Progress: ──────────────────────────────────────────────────────────────────────────────────────────────── 100%
    # Warning: Skipping prob 0.2 because no markers found in both directions
    # Warning: Skipping prob 0.25 because no markers found in both directions
    # Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam,  :
    #   All values in the stats vector are greater than zero and scoreType is "std", maybe you should switch to scoreType = "pos".
    # Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam,  :
    #   All values in the stats vector are greater than zero and scoreType is "std", maybe you should switch to scoreType = "pos".
    # Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam,  :
    #   All values in the stats vector are greater than zero and scoreType is "std", maybe you should switch to scoreType = "pos".
    # Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam,  :
    #   All values in the stats vector are greater than zero and scoreType is "std", maybe you should switch to scoreType = "pos".
    # ✔ [2025/11/07 08:54:19] Optimal threshold: 0.4 (NES difference: 6.442)
    # ℹ [2025/11/07 08:54:19] Valid results: 4/6 thresholds
    # ℹ [2025/11/07 08:54:19] Running fixed threshold mode with prob = 0.4
    # ℹ [2025/11/07 08:54:19] Computing AUC scores...
    # Genes in the gene sets NOT available in the dataset:
    #   gene_pos:   13 (6% of 230)
    #   gene_neg:   54 (12% of 446)
    # ✔ [2025/11/07 08:54:21] Classified 181 Positive, 191 Negative cells
    # ℹ [2025/11/07 08:54:21] Finding markers between `Positive` and `Negative` group...
    # ✔ [2025/11/07 08:54:22] Found 267 positive markers, 2 negative markers
    # ✔ [2025/11/07 08:54:22] scPP screening done.

<!-- ## Find Optimal Parameters for `DEGAS` -->
<!-- ## Find Optimal Parameters for `SCIPAC` -->
<!-- ## Find Optimal Parameters for `SIDISH` -->
<!-- ## Find Optimal Parameters for `TiRank` -->

## Session Info
