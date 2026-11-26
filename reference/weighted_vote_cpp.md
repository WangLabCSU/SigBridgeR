# Combine Votes Using Per-Voter Weights

Combine Votes Using Per-Voter Weights

## Usage

``` r
weighted_vote_cpp(vote_data, weights, ties_method)
```

## Arguments

- vote_data:

  Character matrix of votes, one row per cell and one column per voter.

- weights:

  Numeric vector of weights, one per column of `vote_data`.

- ties_method:

  `integer`. Tie handling: 0 random, 1 `"Positive"`, 2 `"Other"`.

## Value

A character vector of length `nrow(vote_data)` with values `"Positive"`
or `"Other"`.
