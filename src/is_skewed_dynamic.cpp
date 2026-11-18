#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Test Whether a Target Value Is Over- or Under-Represented
//'
//' @param x Vector of values to examine.
//' @param target `numeric`. Value to test. Default: `0`.
//' @param expected_p `numeric`. Expected proportion of `target`.
//'   Default: `0.8`.
//' @param n_sd `integer`. Number of binomial standard errors allowed.
//'   Default: `4L`.
//'
//' @return A length-1 logical: `TRUE` when the observed proportion deviates
//'   by more than `n_sd` standard errors; `NA` for empty input.
//'
//' @keywords internal
// [[Rcpp::export]]
LogicalVector IsSkewedDynamic_cpp(SEXP x, double target = 0.0,
                                  double expected_p = 0.8, int n_sd = 4) {
  R_xlen_t n = Rf_xlength(x);
  if (n == 0) {
    return LogicalVector::create(NA_LOGICAL);
  }
  R_xlen_t count = 0;
  switch (TYPEOF(x)) {
  case REALSXP: {
    double *px = REAL(x);
    for (R_xlen_t i = 0; i < n; ++i) {
      double v = px[i];
      if (Rcpp::NumericVector::is_na(v)) {
        stop("x contains NA");
      }
      if (v == target) {
        ++count;
      }
    }
    break;
  }
  case INTSXP: {
    int *px = INTEGER(x);
    for (R_xlen_t i = 0; i < n; ++i) {
      int v = px[i];
      if (v == NA_INTEGER) {
        stop("x contains NA");
      }
      if (static_cast<double>(v) == target) {
        ++count;
      }
    }
    break;
  }
  case LGLSXP: {
    int *px = LOGICAL(x);
    for (R_xlen_t i = 0; i < n; ++i) {
      int v = px[i];
      if (v == NA_LOGICAL) {
        stop("x contains NA");
      }
      if (static_cast<double>(v) == target) {
        ++count;
      }
    }
    break;
  }
  default:
    stop("x must be numeric, integer, or logical");
  }
  double p_hat = static_cast<double>(count) / static_cast<double>(n);
  double sd_expected =
      std::sqrt(expected_p * (1.0 - expected_p) / static_cast<double>(n));
  bool ans =
      std::abs(p_hat - expected_p) > static_cast<double>(n_sd) * sd_expected;
  return LogicalVector::create(ans);
}
