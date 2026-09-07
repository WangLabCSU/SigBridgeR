#include <RcppArmadillo.h>

#include "SigBridgeR.h"

using namespace Rcpp;

// [[Rcpp::export]]
bool IsCountsMatrixImpl(const arma::mat &x, bool verbose = false,
                        double integer_tol = 1e-8,
                        double min_integer_fraction = 0.95) {
  const arma::uword n = x.n_elem;
  if (n == 0) {
    if (verbose)
      cli_emit("info", "Result: FALSE; matrix is empty");
    return false;
  }
  // Check NA, NaN, and Inf
  if (!x.is_finite()) {
    if (verbose)
      cli_emit("info", "Result: FALSE; matrix contains NA, NaN, or Inf");
    return false;
  }
  // Check non-negativity
  double nonnegative_fraction = static_cast<double>(arma::accu(x >= 0)) / n;
  if (nonnegative_fraction < 1.0) {
    if (verbose)
      cli_emit("info", "Result: FALSE; matrix contains negative values");
    return false;
  }
  // Check whether values are integers
  arma::mat rounded = arma::round(x);
  double integer_fraction =
      static_cast<double>(arma::accu(arma::abs(x - rounded) <= integer_tol)) /
      n;
  if (integer_fraction < min_integer_fraction) {
    if (verbose)
      cli_emit("info", "Result: FALSE; insufficient integer-like values "
                       "(integer fraction " +
                           std::to_string(integer_fraction) + " < threshold " +
                           std::to_string(min_integer_fraction) + ")");
    return false;
  }
  // Check library sizes for all samples
  arma::rowvec library_sizes = arma::sum(x, 0);
  if (arma::any(library_sizes <= 0)) {
    if (verbose)
      cli_emit("info", "Result: FALSE; at least one sample has a library "
                       "size <= 0");
    return false;
  }
  if (verbose) {
    double zero_fraction = static_cast<double>(arma::accu(x == 0)) / n;
    cli_emit("success",
             "Result: TRUE; matrix is likely a counts matrix "
             "(non-negative fraction " +
                 std::to_string(nonnegative_fraction) + ", integer fraction " +
                 std::to_string(integer_fraction) + ", zero fraction " +
                 std::to_string(zero_fraction) + ")");
  }
  return true;
}
