#' Test Whether the Proportion of a Target Value Is Highly Skewed
#'
#' @description
#' `IsSkewedDynamic()` tests whether the observed proportion of a specified
#' target/background value in a vector deviates substantially from an expected
#' proportion.
#'
#' The function is designed for detecting highly skewed phenotype-like vectors,
#' especially binary or discrete vectors where one value represents the
#' background/control class. If the observed proportion of `target` differs from
#' `expected_p` by more than `n_sd` standard deviations under a binomial model,
#' the vector is considered skewed.
#'
#' Formally, let:
#'
#' \deqn{\hat{p} = \frac{\#(x_i = target)}{n}}
#'
#' and
#'
#' \deqn{SD = \sqrt{\frac{expected\_p(1 - expected\_p)}{n}}}
#'
#' The function returns `TRUE` when:
#'
#' \deqn{|\hat{p} - expected\_p| > n\_sd \times SD}
#'
#' Otherwise, it returns `FALSE`.
#'
#' If `x` has length zero, the function returns `NA`.
#'
#' @param x A vector to be tested. Must be numeric, integer, or logical.
#'   Missing values are not allowed.
#' @param target A single numeric value indicating the target/background value
#'   whose proportion should be evaluated. Default is `0`.
#' @param expected_p A single numeric value between 0 and 1 indicating the
#'   expected proportion of `target` in `x`. Default is `0.8`.
#' @param n_sd A single integer or numeric value specifying how many standard
#'   deviations away from `expected_p` are allowed before the vector is
#'   considered skewed. Default is `4L`.
#'
#' @return
#' A logical vector of length 1:
#'
#' * `TRUE` if the observed proportion of `target` deviates from `expected_p`
#'   by more than `n_sd` standard deviations.
#' * `FALSE` otherwise.
#' * `NA` if `x` has length zero.
#'
#' @details
#' This function assumes the expected occurrence of `target` follows a binomial
#' proportion model with probability `expected_p`. The standard deviation is
#' computed as:
#'
#' \deqn{\sqrt{expected\_p(1 - expected\_p) / n}}
#'
#' where `n` is the length of `x`.
#'
#' The check is symmetric: both an unexpectedly high and unexpectedly low
#' proportion of `target` are considered skewed.
#'
#' Internally, the computation is delegated to a C++ implementation via
#' Rcpp for efficiency.
#'
#' @section Missing values:
#' `x` must not contain missing values. If any `NA` is present, the function
#' throws an error.
#'
#' @section Supported input types:
#' `x` may be one of:
#'
#' * numeric
#' * integer
#' * logical
#'
#' Other input types are rejected by the underlying C++ implementation.
#'
#' @examples
#' x1 <- c(rep(0, 80), rep(1, 20))
#' IsSkewedDynamic(x1)
#'
#' x2 <- c(rep(0, 95), rep(1, 5))
#' IsSkewedDynamic(x2)
#'
#' x3 <- c(rep(0, 50), rep(1, 50))
#' IsSkewedDynamic(x3, expected_p = 0.8)
#'
#' IsSkewedDynamic(logical(c(rep(FALSE, 80), rep(TRUE, 20))),
#'                   target = 0,
#'                   expected_p = 0.8)
#'
#' IsSkewedDynamic(numeric(0))
#'
#' @export
IsSkewedDynamic <- function(
  x,
  target = 0,
  expected_p = 0.8,
  n_sd = 4L
) {
  if (anyNA(x)) {
    Abort("x contains {.val NA}")
  }
  IsSkewedDynamic_cpp(x, target, expected_p, n_sd)
}
