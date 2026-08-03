# WeightedVote - basic voting --------------------------------------------

describe("WeightedVote - basic voting", {
  # arg_match() is called without rlang:: prefix in the source; mock it.
  local_mocked_bindings(
    arg_match = rlang::arg_match,
    .package = "SigBridgeR"
  )

  it("returns Positive when all voters vote Positive", {
    vote_data <- data.frame(
      V1 = "Positive",
      V2 = "Positive",
      V3 = "Positive"
    )
    weights <- c(V1 = 1, V2 = 1, V3 = 1)

    result <- WeightedVote(vote_data, weights)

    expect_equal(unname(result), "Positive")
  })

  it("returns Other when all voters vote Negative/Neutral/Other", {
    vote_data <- data.frame(
      V1 = "Negative",
      V2 = "Neutral",
      V3 = "Other"
    )
    weights <- c(V1 = 1, V2 = 1, V3 = 1)

    result <- WeightedVote(vote_data, weights)

    expect_equal(unname(result), "Other")
  })

  it("majority Positive wins with equal weights", {
    vote_data <- data.frame(
      V1 = c("Positive", "Neutral", "Positive"),
      V2 = c("Positive", "Neutral", "Negative"),
      V3 = c("Negative", "Positive", "Positive")
    )
    weights <- c(V1 = 1, V2 = 1, V3 = 1)

    result <- WeightedVote(vote_data, weights)

    # Row 1: 2 Positive, 1 Negative -> Positive
    # Row 2: 1 Positive (V3), 2 Other (Neutral x2) -> Other
    # Row 3: 2 Positive, 1 Negative -> Positive
    expect_equal(unname(result), c("Positive", "Other", "Positive"))
  })

  it("majority Other wins with equal weights", {
    vote_data <- data.frame(
      V1 = c("Negative", "Positive", "Negative"),
      V2 = c("Neutral", "Negative", "Other"),
      V3 = c("Other", "Neutral", "Neutral")
    )
    weights <- c(V1 = 1, V2 = 1, V3 = 1)

    result <- WeightedVote(vote_data, weights)

    # Row 1: 1 Negative + 1 Neutral + 1 Other -> 3 Other -> Other
    # Row 2: 1 Positive vs 1 Negative + 1 Neutral -> 2 Other -> Other
    # Row 3: 1 Negative + 1 Other + 1 Neutral -> 3 Other -> Other
    expect_equal(unname(result), c("Other", "Other", "Other"))
  })

  it("single voter works", {
    vote_data <- data.frame(
      Expert = c("Positive", "Negative", "Neutral", "Other")
    )
    weights <- c(Expert = 1)

    result <- WeightedVote(vote_data, weights)

    expect_equal(unname(result), c("Positive", "Other", "Other", "Other"))
  })

  it("single row with multiple voters", {
    vote_data <- data.frame(
      V1 = "Positive",
      V2 = "Negative",
      V3 = "Neutral"
    )
    weights <- c(V1 = 1, V2 = 1, V3 = 1)

    result <- WeightedVote(vote_data, weights)

    # 1 Positive vs 2 Other -> Other wins
    expect_equal(unname(result), "Other")
  })

  it("returns character vector", {
    vote_data <- data.frame(
      V1 = c("Positive", "Negative"),
      V2 = c("Positive", "Positive")
    )
    weights <- c(V1 = 1, V2 = 1)

    result <- WeightedVote(vote_data, weights)

    expect_type(result, "character")
    expect_equal(length(result), 2)
  })
})

# WeightedVote - weighted scoring ----------------------------------------

describe("WeightedVote - weighted scoring", {
  local_mocked_bindings(
    arg_match = rlang::arg_match,
    .package = "SigBridgeR"
  )

  it("higher weight overrides majority", {
    vote_data <- data.frame(
      V1 = c("Positive", "Negative"),
      V2 = c("Negative", "Positive"),
      V3 = c("Negative", "Positive")
    )
    # V1 has much higher weight
    weights <- c(V1 = 100, V2 = 1, V3 = 1)

    result <- WeightedVote(vote_data, weights)

    # Row 1: V1 Positive (100) vs V2+V3 Other (2) -> Positive
    # Row 2: V1 Negative (100) vs V2+V3 Positive (2) -> Other
    expect_equal(unname(result), c("Positive", "Other"))
  })

  it("non-integer weights work correctly", {
    vote_data <- data.frame(
      V1 = c("Positive", "Negative"),
      V2 = c("Negative", "Negative")
    )
    weights <- c(V1 = 1.5, V2 = 1.0)

    result <- WeightedVote(vote_data, weights)

    # Row 1: Positive 1.5 > Other 1.0 -> Positive
    # Row 2: Other 1.5 > Other 1.0 -> Other (no Positive votes)
    expect_equal(unname(result), c("Positive", "Other"))
  })

  it("zero weight voter has no influence", {
    vote_data <- data.frame(
      V1 = c("Positive", "Negative"),
      V2 = c("Negative", "Negative")
    )
    weights <- c(V1 = 0, V2 = 1)

    result <- WeightedVote(vote_data, weights)

    # Row 1: V1 has weight 0, V2 says Negative -> Other
    # Row 2: both Negative -> Other
    expect_equal(unname(result), c("Other", "Other"))
  })

  it("weights in different order than vote_data columns are reordered", {
    vote_data <- data.frame(
      Expert1 = c("Positive", "Negative"),
      Expert2 = c("Negative", "Negative")
    )
    # weights in reverse column order
    weights <- c(Expert2 = 100, Expert1 = 1)

    result <- WeightedVote(vote_data, weights)

    # Expert2 has higher weight, says Negative for both -> Other
    expect_equal(unname(result), c("Other", "Other"))
  })

  it("very small weights still contribute", {
    vote_data <- data.frame(
      V1 = c("Positive"),
      V2 = c("Negative")
    )
    weights <- c(V1 = 1e-10, V2 = 1e-11)

    result <- WeightedVote(vote_data, weights)

    expect_equal(unname(result), "Positive")
  })
})

# WeightedVote - ties.method ---------------------------------------------

describe("WeightedVote - ties.method", {
  local_mocked_bindings(
    arg_match = rlang::arg_match,
    .package = "SigBridgeR"
  )

  it("ties.method = 'first' chooses Positive on tie", {
    # TieFirst::choose() returns 0 -> "Positive"
    vote_data <- data.frame(
      V1 = c("Positive", "Negative"),
      V2 = c("Negative", "Positive")
    )
    weights <- c(V1 = 1, V2 = 1)

    result <- WeightedVote(vote_data, weights, ties.method = "first")

    # Both rows are ties: first wins -> always Positive
    expect_equal(unname(result), c("Positive", "Positive"))
  })

  it("ties.method = 'last' chooses Other on tie", {
    # TieLast::choose() returns 1 -> "Other"
    vote_data <- data.frame(
      V1 = c("Positive", "Negative"),
      V2 = c("Negative", "Positive")
    )
    weights <- c(V1 = 1, V2 = 1)

    result <- WeightedVote(vote_data, weights, ties.method = "last")

    # Both rows are ties: last wins -> always Other
    expect_equal(unname(result), c("Other", "Other"))
  })

  it("ties.method = 'random' produces valid output", {
    vote_data <- data.frame(
      V1 = c("Positive", "Negative"),
      V2 = c("Negative", "Positive")
    )
    weights <- c(V1 = 1, V2 = 1)

    result <- WeightedVote(vote_data, weights, ties.method = "random")

    # Output should be either Positive or Other
    expect_true(all(result %in% c("Positive", "Other")))
    expect_equal(length(result), 2)
  })

  it("ties.method is case-sensitive and rejects wrong case", {
    vote_data <- data.frame(V1 = c("Positive"), V2 = c("Negative"))
    weights <- c(V1 = 1, V2 = 1)

    expect_error(
      WeightedVote(vote_data, weights, ties.method = "Random"),
      class = "rlang_error"
    )
  })

  it("invalid ties.method aborts", {
    vote_data <- data.frame(V1 = c("Positive"), V2 = c("Negative"))
    weights <- c(V1 = 1, V2 = 1)

    expect_error(
      WeightedVote(vote_data, weights, ties.method = "unknown"),
      class = "rlang_error"
    )
  })
})

# WeightedVote - row names -----------------------------------------------

describe("WeightedVote - row names", {
  local_mocked_bindings(
    arg_match = rlang::arg_match,
    .package = "SigBridgeR"
  )

  it("preserves row names in output", {
    vote_data <- data.frame(
      Expert1 = c("Positive", "Negative", "Neutral"),
      Expert2 = c("Positive", "Neutral", "Other"),
      row.names = c("Cell1", "Cell2", "Cell3")
    )
    weights <- c(Expert1 = 1, Expert2 = 1)

    result <- WeightedVote(vote_data, weights)

    expect_equal(names(result), c("Cell1", "Cell2", "Cell3"))
    expect_equal(
      unname(result),
      c("Positive", "Other", "Other")
    )
  })

  it("returns unnamed vector when no row names", {
    # Use matrix without dimnames to avoid auto row names from data.frame
    vote_mat <- matrix(
      c("Positive", "Negative", "Positive", "Positive"),
      nrow = 2,
      ncol = 2
    )
    colnames(vote_mat) <- c("V1", "V2")
    weights <- c(V1 = 1, V2 = 1)

    result <- WeightedVote(vote_mat, weights)

    expect_null(names(result))
  })

  it("matrix input with row names preserves them", {
    vote_mat <- matrix(
      c("Positive", "Negative", "Neutral", "Other"),
      nrow = 2,
      ncol = 2,
      dimnames = list(c("R1", "R2"), c("V1", "V2"))
    )
    weights <- c(V1 = 1, V2 = 1)

    result <- WeightedVote(vote_mat, weights)

    expect_equal(names(result), c("R1", "R2"))
  })
})

# WeightedVote - edge cases ----------------------------------------------

describe("WeightedVote - edge cases", {
  local_mocked_bindings(
    arg_match = rlang::arg_match,
    .package = "SigBridgeR"
  )

  it("handles data with only Negative/Neutral/Other (no Positive)", {
    vote_data <- data.frame(
      V1 = c("Negative", "Neutral"),
      V2 = c("Other", "Negative")
    )
    weights <- c(V1 = 1, V2 = 1)

    result <- WeightedVote(vote_data, weights)

    expect_equal(unname(result), c("Other", "Other"))
  })

  it("handles data with only Positive", {
    vote_data <- data.frame(
      V1 = c("Positive", "Positive"),
      V2 = c("Positive", "Positive")
    )
    weights <- c(V1 = 1, V2 = 1)

    result <- WeightedVote(vote_data, weights)

    expect_equal(unname(result), c("Positive", "Positive"))
  })

  it("matrix input works the same as data.frame", {
    df <- data.frame(
      V1 = c("Positive", "Negative"),
      V2 = c("Negative", "Positive"),
      V3 = c("Positive", "Positive")
    )
    mat <- as.matrix(df)
    weights <- c(V1 = 1, V2 = 1, V3 = 1)

    result_df <- WeightedVote(df, weights)
    result_mat <- WeightedVote(mat, weights)

    expect_equal(unname(result_df), unname(result_mat))
  })

  it("NA votes are ignored (do not contribute to either side)", {
    vote_data <- data.frame(
      V1 = c("Positive", NA),
      V2 = c("Negative", "Positive")
    )
    weights <- c(V1 = 1, V2 = 1)

    # Use first to avoid non-deterministic random tie-breaking
    result <- WeightedVote(vote_data, weights, ties.method = "first")

    # Row 1: V1 Positive (1), V2 Other (1) -> tie, first wins -> Positive
    # Row 2: V1 NA (ignored), V2 Positive (1) -> Positive
    expect_equal(unname(result), c("Positive", "Positive"))
  })

  it("NaN weights are treated as zero (no contribution)", {
    vote_data <- data.frame(
      V1 = c("Positive"),
      V2 = c("Negative")
    )
    weights <- c(V1 = NaN, V2 = 1)

    result <- WeightedVote(vote_data, weights)

    # V1 has NaN weight (ignored by ISNAN check), V2 says Negative -> Other
    expect_equal(unname(result), "Other")
  })

  it("all-zero weights with ties.method = 'first' returns Positive", {
    vote_data <- data.frame(
      V1 = c("Positive"),
      V2 = c("Negative")
    )
    weights <- c(V1 = 0, V2 = 0)

    result <- WeightedVote(vote_data, weights, ties.method = "first")

    # 0 vs 0 tie, TieFirst::choose() returns 0 -> Positive
    expect_equal(unname(result), "Positive")
  })

  it("all-zero weights with ties.method = 'last' returns Other", {
    vote_data <- data.frame(
      V1 = c("Positive"),
      V2 = c("Negative")
    )
    weights <- c(V1 = 0, V2 = 0)

    result <- WeightedVote(vote_data, weights, ties.method = "last")

    # 0 vs 0 tie, TieLast::choose() returns 1 -> Other
    expect_equal(unname(result), "Other")
  })

  it("many rows and voters works", {
    n_rows <- 50
    n_voters <- 10
    set.seed(42)

    labels <- c("Positive", "Negative", "Neutral", "Other")
    vote_data <- as.data.frame(
      matrix(
        sample(labels, n_rows * n_voters, replace = TRUE),
        nrow = n_rows,
        ncol = n_voters
      )
    )
    colnames(vote_data) <- paste0("V", seq_len(n_voters))
    weights <- setNames(rep(1, n_voters), paste0("V", seq_len(n_voters)))

    result <- WeightedVote(vote_data, weights)

    expect_equal(length(result), n_rows)
    expect_true(all(result %in% c("Positive", "Other")))
  })
})

# WeightedVote - input validation ----------------------------------------

describe("WeightedVote - input validation", {
  local_mocked_bindings(
    arg_match = rlang::arg_match,
    .package = "SigBridgeR"
  )

  it("aborts when vote_data has no column names", {
    vote_data <- matrix(c("Positive", "Negative"), nrow = 1)
    # matrix without dimnames
    weights <- c(V1 = 1, V2 = 1)

    expect_error(
      WeightedVote(vote_data, weights),
      "must have column names"
    )
  })

  it("aborts when weights names do not match vote_data columns", {
    vote_data <- data.frame(
      Expert1 = c("Positive", "Negative"),
      Expert2 = c("Negative", "Positive")
    )
    weights <- c(Wrong1 = 1, Wrong2 = 1)

    expect_error(
      WeightedVote(vote_data, weights),
      "must have the same names"
    )
  })

  it("aborts when weights have different length than columns", {
    vote_data <- data.frame(
      V1 = c("Positive", "Negative"),
      V2 = c("Negative", "Positive"),
      V3 = c("Positive", "Positive")
    )
    weights <- c(V1 = 1, V2 = 1)

    expect_error(
      WeightedVote(vote_data, weights),
      "must have the same names"
    )
  })

  it("aborts when weights have duplicate names", {
    vote_data <- data.frame(
      V1 = c("Positive"),
      V2 = c("Negative")
    )
    weights <- c(V1 = 1, V1 = 2)

    expect_error(
      WeightedVote(vote_data, weights),
      "must have the same names"
    )
  })

  it("aborts when weights are not numeric", {
    vote_data <- data.frame(V1 = c("Positive"), V2 = c("Negative"))
    weights <- c(V1 = "high", V2 = "low")

    expect_error(
      WeightedVote(vote_data, weights),
      class = "chk_error"
    )
  })
})
