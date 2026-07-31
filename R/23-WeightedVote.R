#' @title Weighted Voting Aggregation for Multi-Voter Classification
#' @description
#' Aggregates classification votes from multiple voters using weighted scoring.
#' Each voter assigns one of four labels ("Positive", "Negative", "Neutral", "Other")
#' to each item, and votes are combined using user-specified weights.
#'
#' @param vote_data A data frame or matrix where:
#'   \itemize{
#'     \item Rows represent items to be classified
#'     \item Columns represent individual voters
#'     \item Cell values must be one of: "Positive", "Negative", "Neutral", "Other"
#'   }
#'   Row names (if present) are preserved in the output.
#' @param weights Named numeric vector specifying the weight for each voter.
#'   Names must exactly match \code{colnames(vote_data)}. Weights should be
#'   non-negative, with higher values indicating greater influence.
#' @param ties.method Character string specifying how to break ties.
#'   Must be one of "random", "first", or "last".
#'
#' @return Character vector of final aggregated labels, one per row in \code{vote_data}.
#'   Each element is one of "Positive" or "Other" (the winner with
#'   highest weighted score). If input has row names, they are preserved as names
#'   in the output vector.
#'
#' @export
#' @examples
#' \donttest{
#' # Example 1: Simple majority voting
#' vote_data <- data.frame(
#'   Expert1 = c("Positive", "Negative", "Neutral", "Positive"),
#'   Expert2 = c("Positive", "Neutral", "Other", "Negative"),
#'   Expert3 = c("Negative", "Positive", "Neutral", "Positive"),
#'   row.names = c("Cell1", "Cell2", "Cell3", "Cell4")
#' )
#'
#' weights <- c(Expert1 = 1, Expert2 = 1, Expert3 = 1)
#' result <- WeightedVote(vote_data, weights)
#' }
#'
#'
#'
WeightedVote <- function(
  vote_data = data.frame(),
  weights = double(),
  ties.method = c("random", "first", "last")
) {
  ties.method <- arg_match(ties.method)
  voter_cols <- colnames(vote_data)
  weight_names <- names(weights)

  chk::chk_numeric(weights)

  if (is.null(voter_cols)) {
    Abort("`vote_data` must have column names.")
  }

  if (
    length(weights) != length(voter_cols) ||
      anyDuplicated(weight_names) ||
      anyDuplicated(voter_cols) ||
      !setequal(weight_names, voter_cols)
  ) {
    Abort(
      "The `weights` must have the same names as `vote_data`",
      "Colnames of `vote_data`: {.val {voter_cols}}",
      "Colnames of `weights`: {.val {weight_names}}"
    )
  }

  w <- unname(weights[voter_cols])

  vote_mat <- as.matrix(vote_data)
  storage.mode(vote_mat) <- "character"

  ties_code <- switch(
    ties.method,
    random = 0L,
    first = 1L,
    last = 2L
  )

  res <- weighted_vote_cpp(
    vote_data = vote_mat,
    weights = w,
    ties_method = ties_code
  )

  r_n <- rownames(vote_data)
  if (!is.null(r_n)) {
    names(res) <- r_n
  }

  res
}
