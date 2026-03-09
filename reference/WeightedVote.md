# Weighted Voting Aggregation for Multi-Voter Classification

Aggregates classification votes from multiple voters using weighted
scoring. Each voter assigns one of four labels ("Positive", "Negative",
"Neutral", "Other") to each item, and votes are combined using
user-specified weights. The "Other" category contributes half-weight to
both "Negative" and "Neutral" categories.

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
Each element is one of "Positive", "Negative", or "Neutral" (the winner
with highest weighted score). If input has row names, they are preserved
as names in the output vector.

## Scoring Logic

For each item (row), scores are computed as:

- Positive score:

  Sum of weights where vote = "Positive"

- Negative score:

  Sum of weights where vote = "Negative" + 0.5 × sum of weights where
  vote = "Other"

- Neutral score:

  Sum of weights where vote = "Neutral" + 0.5 × sum of weights where
  vote = "Other"

The category with the highest score wins. In case of ties, the first
category (in order: Positive, Negative, Neutral) is selected.

## See also

[`table`](https://rdrr.io/r/base/table.html) for vote counting,
[`max.col`](https://rdrr.io/r/base/maxCol.html) for winner selection

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Simple majority voting
vote_data <- data.frame(
  Expert1 = c("Positive", "Negative", "Neutral", "Positive"),
  Expert2 = c("Positive", "Neutral", "Other", "Negative"),
  Expert3 = c("Negative", "Positive", "Neutral", "Positive"),
  row.names = c("Gene1", "Gene2", "Gene3", "Gene4")
)

weights <- c(Expert1 = 1, Expert2 = 1, Expert3 = 1)
result <- WeightedVote(vote_data, weights)
# Gene1: Positive(2) vs Negative(1) → "Positive"
# Gene2: each has 1 vote → "Positive" (first in tie)
# Gene3: Neutral(1.5) vs Other(0.5) → "Neutral"
# Gene4: Positive(2) vs Negative(1) → "Positive"

# Example 2: Unequal weights (expert confidence levels)
weights <- c(Expert1 = 3, Expert2 = 1, Expert3 = 2)
result <- WeightedVote(vote_data, weights)

# Example 3: "Other" category splits weight
vote_data <- data.frame(
  Voter1 = c("Positive", "Other", "Negative"),
  Voter2 = c("Other", "Neutral", "Other")
)
weights <- c(Voter1 = 1, Voter2 = 1)
result <- WeightedVote(vote_data, weights)
# Row1: Positive(1) + Other(0.5 to Neg/Neu) → Positive wins
# Row2: Other(0.5 to Neg/Neu) + Neutral(1) → Neutral wins
# Row3: Negative(1) + Other(0.5 to Neg/Neu) → Negative wins

# Example 4: Programmatic use with row names
result <- WeightedVote(vote_data, weights)
names(result)  # Returns row names from vote_data
table(result)   # Count votes per category
} # }
```
