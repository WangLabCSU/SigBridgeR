#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// Min-max normalization.
// - NA / NaN / +-Inf 视为缺失：不参与 min/max，输出 NA
// - 有效值全部相同（或仅 1 个有效值）：有效位置输出 0.5
// [[Rcpp::export]]
NumericVector normalize_metric(const NumericVector &x, bool invert = false) {
  const R_xlen_t n = x.size();
  const double *px = x.begin();

  double xmin = R_PosInf;
  double xmax = R_NegInf;
  for (R_xlen_t i = 0; i < n; ++i) {
    const double xi = px[i];
    if (!std::isfinite(xi))
      continue;
    if (xi < xmin)
      xmin = xi;
    if (xi > xmax)
      xmax = xi;
  }

  NumericVector out = no_init(n);
  double *po = out.begin();

  // 无有效值 / 仅一个有效值 / 全部相同
  if (!(xmax > xmin)) {
    for (R_xlen_t i = 0; i < n; ++i) {
      po[i] = std::isfinite(px[i]) ? 0.5 : NA_REAL;
    }
    return out;
  }

  const double range = xmax - xmin;
  if (invert) {
    for (R_xlen_t i = 0; i < n; ++i) {
      const double xi = px[i];
      po[i] = std::isfinite(xi) ? (xmax - xi) / range : NA_REAL;
    }
  } else {
    for (R_xlen_t i = 0; i < n; ++i) {
      const double xi = px[i];
      po[i] = std::isfinite(xi) ? (xi - xmin) / range : NA_REAL;
    }
  }
  return out;
}