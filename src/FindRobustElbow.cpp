#include "Rtatami.h"
#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional> // std::greater
#include <vector>

// [[Rcpp::depends(beachmat, assorthead)]]
// [[Rcpp::plugins(cpp17)]]

using namespace Rcpp;

/*
 * Numeric core of FindRobustElbow().
 *
 * initialized_matrix must be the external pointer returned on the R side by:
 *
 *     beachmat::initializeCpp(obj[["pca"]]@cell.embeddings)
 *
 * i.e. a cells x PCs matrix whose columns are the PCA embeddings. The columns
 * are streamed through tatami, so none of the R-side matrix classes need to be
 * copied or coerced to a dense base matrix first.
 *
 * The first `ndims` columns are summarised as:
 *
 *   stdev               = column standard deviation
 *   variance            = stdev^2
 *   pct_variance        = variance / sum(variance) * 100
 *   cumulative_variance = cumsum(pct_variance)
 *
 * and then fed into the same six heuristics as the original R implementation.
 *
 * Package notes:
 *   DESCRIPTION: LinkingTo: Rcpp, beachmat, assorthead
 *                SystemRequirements: C++17
 *   src/Makevars(.win): CXX_STD = CXX17
 */
namespace {

// Default cut-off reported by the "last index above threshold" heuristics when
// no PC qualifies. Always clamped to the number of PCs actually used.
constexpr int kDefaultFallback = 10;

// Relative tolerance for the strict ">" comparisons against the mean / 2*sd
// thresholds. Without it, near-equal PC variances make the outcome depend on
// rounding noise (e.g. every pct_variance == mean up to 1 ulp).
constexpr double kRelTol = 1e-12;

// Sample standard deviation (n - 1 denominator) of the first `ndims` columns.
// Precondition (checked): 1 <= ndims <= ncol and nrow >= 2.
std::vector<double> column_sd(const tatami::NumericMatrix *matrix, int ndims) {
  const int nr = matrix->nrow();
  const int nc = matrix->ncol();
  if (ndims < 1 || ndims > nc) {
    Rcpp::stop("Internal error: 'ndims' (%d) is outside [1, %d].", ndims, nc);
  }
  if (nr < 2) {
    Rcpp::stop("Internal error: at least 2 rows are required.");
  }

  std::vector<double> sdev(static_cast<std::size_t>(ndims), 0.0);
  std::vector<double> buffer(static_cast<std::size_t>(nr));

  auto ext = matrix->dense_column();
  const double n = static_cast<double>(nr);
  const double n_minus_1 = static_cast<double>(nr - 1);

  for (int c = 0; c < ndims; ++c) {
    Rcpp::checkUserInterrupt();
    const double *values = ext->fetch(c, buffer.data());

    // Pass 1: mean. NaN / +-Inf propagate into `sum`, so a single check after
    // the loop replaces the per-element branch.
    double sum = 0.0;
    for (int r = 0; r < nr; ++r)
      sum += values[r];
    if (!std::isfinite(sum)) {
      Rcpp::stop("PCA embeddings contain non-finite values (PC %d).", c + 1);
    }
    const double mean = sum / n;

    // Pass 2: sum of squared deviations (two-pass for numerical stability).
    double ss = 0.0;
    for (int r = 0; r < nr; ++r) {
      const double d = values[r] - mean;
      ss += d * d;
    }

    sdev[static_cast<std::size_t>(c)] = std::sqrt(ss / n_minus_1);
  }

  return sdev;
}

double mean_of(const std::vector<double> &x) {
  if (x.empty())
    return NA_REAL;
  double sum = 0.0;
  for (const double v : x)
    sum += v;
  return sum / static_cast<double>(x.size());
}

// Sample standard deviation (n - 1 denominator), matching stats::sd().
double sd_of(const std::vector<double> &x) {
  if (x.size() < 2)
    return NA_REAL;
  const double mean = mean_of(x);
  double ss = 0.0;
  for (const double v : x) {
    const double d = v - mean;
    ss += d * d;
  }
  return std::sqrt(ss / static_cast<double>(x.size() - 1));
}

// 1-based index of the first element strictly above `threshold`.
int first_index_above(const std::vector<double> &x, double threshold,
                      int fallback) {
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (x[i] > threshold)
      return static_cast<int>(i) + 1;
  }
  return fallback;
}

// 1-based index of the last element strictly above `threshold`, where the
// threshold is nudged upwards by a relative tolerance so that ties up to
// rounding error are not counted as "above". A NaN threshold never matches
// and yields `fallback`.
int last_index_above(const std::vector<double> &x, double threshold,
                     int fallback) {
  const double cut = threshold + kRelTol * std::fabs(threshold);
  for (std::size_t i = x.size(); i-- > 0;) {
    if (x[i] > cut)
      return static_cast<int>(i) + 1;
  }
  return fallback;
}

// Elbow from the largest absolute second derivative of the variance curve.
int second_derivative_elbow(const std::vector<double> &variance, int ndims) {
  if (ndims < 3)
    return ndims;

  std::vector<double> first(static_cast<std::size_t>(ndims - 1), 0.0);
  for (int i = 0; i < ndims - 1; ++i)
    first[static_cast<std::size_t>(i)] =
        variance[static_cast<std::size_t>(i) + 1] -
        variance[static_cast<std::size_t>(i)];

  std::size_t best = 0;
  double best_value = std::fabs(first[1] - first[0]);
  for (std::size_t k = 1; k + 1 < first.size(); ++k) {
    const double value = std::fabs(first[k + 1] - first[k]);
    if (value > best_value) {
      best_value = value;
      best = k;
    }
  }

  return std::min(static_cast<int>(best) + 2, ndims);
}

// Elbow as the point farthest from the chord spanning the first and last PC.
int distance_elbow(const std::vector<double> &variance, int ndims) {
  if (ndims < 2)
    return 1;

  const double y_start = variance.front();
  const double y_end = variance[static_cast<std::size_t>(ndims - 1)];
  const double dy = y_end - y_start;
  const double dx = static_cast<double>(ndims - 1);
  const double denominator = std::sqrt(dy * dy + dx * dx);
  if (!(denominator > 0.0))
    return 1;

  int best = 1;
  double best_value = -1.0;
  for (int i = 0; i < ndims; ++i) {
    const double x = static_cast<double>(i + 1);
    const double y = variance[static_cast<std::size_t>(i)];
    const double numerator = std::fabs(
        dy * x - dx * y + static_cast<double>(ndims) * y_start - y_end);
    const double value = numerator / denominator;
    if (value > best_value) {
      best_value = value;
      best = i + 1;
    }
  }

  return best;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List find_robust_elbow_cpp(SEXP initialized_matrix, int ndims) {
  if (ndims < 1)
    Rcpp::stop("'ndims' must be a positive integer.");

  Rtatami::BoundNumericPointer parsed(initialized_matrix);
  const auto &matrix = parsed->ptr;

  const int nc = matrix->ncol();
  if (nc < 1)
    Rcpp::stop("The PCA embeddings matrix has no columns.");
  if (matrix->nrow() < 2) {
    Rcpp::stop("At least 2 cells are required to compute the standard "
               "deviation of each PC.");
  }

  // The R wrapper already warns about and clamps this; keep the guard so the
  // kernel is safe when called directly.
  ndims = std::min(ndims, nc);

  const std::vector<double> sdev = column_sd(matrix.get(), ndims);

  const std::size_t n = static_cast<std::size_t>(ndims);
  std::vector<double> variance(n, 0.0);
  double total_variance = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    variance[i] = sdev[i] * sdev[i];
    total_variance += variance[i];
  }
  if (!(total_variance > 0.0) || !std::isfinite(total_variance)) {
    Rcpp::stop("The selected PCs have zero or non-finite total variance; "
               "cannot select a number of dimensions.");
  }

  std::vector<double> pct_variance(n, 0.0);
  std::vector<double> cumulative_variance(n, 0.0);
  double running = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    pct_variance[i] = variance[i] / total_variance * 100.0;
    running += pct_variance[i];
    cumulative_variance[i] = running;
  }

  // Fallback must never exceed the number of PCs actually used, otherwise the
  // reported cut-off (and the plot's vertical line) lies outside the data.
  const int fallback = std::min(kDefaultFallback, ndims);

  // Method 1: variance-based heuristics.
  const int cumulative_90 = first_index_above(cumulative_variance, 90.0, ndims);
  const int cumulative_80 = first_index_above(cumulative_variance, 80.0, ndims);
  const int above_mean =
      last_index_above(pct_variance, mean_of(pct_variance), fallback);
  const int above_2sd =
      last_index_above(pct_variance, 2.0 * sd_of(pct_variance), fallback);

  // Method 2: second derivative.
  const int second_deriv = second_derivative_elbow(variance, ndims);

  // Method 3: distance to the first-last chord.
  const int distance = distance_elbow(variance, ndims);

  // The second largest candidate is the recommended number of dimensions.
  std::vector<int> candidates = {cumulative_90,
                                 cumulative_80,
                                 above_mean,
                                 above_2sd,
                                 second_deriv,
                                 distance,
                                 10,
                                 15,
                                 20,
                                 25,
                                 30,
                                 35};
  std::sort(candidates.begin(), candidates.end(), std::greater<int>());
  // The fixed anchors can exceed the number of PCs that actually exist.
  const int final_dims = std::min(candidates[1], ndims);

  return Rcpp::List::create(
      Rcpp::Named("stdev") = sdev, Rcpp::Named("variance") = variance,
      Rcpp::Named("pct_variance") = pct_variance,
      Rcpp::Named("cumulative_variance") = cumulative_variance,
      Rcpp::Named("cumulative_90") = cumulative_90,
      Rcpp::Named("cumulative_80") = cumulative_80,
      Rcpp::Named("above_mean") = above_mean,
      Rcpp::Named("above_2sd") = above_2sd,
      Rcpp::Named("second_deriv") = second_deriv,
      Rcpp::Named("distance") = distance,
      Rcpp::Named("final_dims") = final_dims);
}