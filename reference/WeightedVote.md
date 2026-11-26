# Weighted Voting Aggregation for Multi-Voter Classification

Aggregates classification votes from multiple voters using weighted
scoring. Each voter assigns one of four labels ("Positive", "Negative",
"Neutral", "Other") to each item, and votes are combined using
user-specified weights.

## Usage

``` r
WeightedVote(
  vote_data = data.frame(),
  weights = double(),
  ties.method = c("random", "first", "last")
)
```

## Arguments

- vote_data:

  A data frame or matrix where:

  - Rows represent items to be classified

  - Columns represent individual voters

  - Cell values must be one of: "Positive", "Negative", "Neutral",
    "Other"

  Row names (if present) are preserved in the output.

- weights:

  Named numeric vector specifying the weight for each voter. Names must
  exactly match `colnames(vote_data)`. Weights should be non-negative,
  with higher values indicating greater influence.

- ties.method:

  Character string specifying how to break ties. Must be one of
  "random", "first", or "last".

## Value

Character vector of final aggregated labels, one per row in `vote_data`.
Each element is one of "Positive" or "Other" (the winner with highest
weighted score). If input has row names, they are preserved as names in
the output vector.

## Examples

``` r
# \donttest{
# Example 1: Simple majority voting
vote_data <- data.frame(
  Expert1 = c("Positive", "Negative", "Neutral", "Positive"),
  Expert2 = c("Positive", "Neutral", "Other", "Negative"),
  Expert3 = c("Negative", "Positive", "Neutral", "Positive"),
  row.names = c("Cell1", "Cell2", "Cell3", "Cell4")
)

weights <- c(Expert1 = 1, Expert2 = 1, Expert3 = 1)
result <- WeightedVote(vote_data, weights)
# }


```
