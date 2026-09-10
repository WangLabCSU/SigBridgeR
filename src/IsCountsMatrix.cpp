#include "Rtatami.h"
#include <Rcpp.h>

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "SigBridgeR.h"

// [[Rcpp::depends(beachmat, assorthead)]]
// [[Rcpp::plugins(cpp17)]]

using namespace Rcpp;

/*
 * Determine whether an expression matrix is likely an RNA-seq counts matrix.
 *
 * initialized_matrix must be the external pointer returned on the R side by:
 *
 *     beachmat::initializeCpp(x)
 *
 * Checks performed:
 *
 * 1. the matrix is non-empty;
 * 2. it contains no NA, NaN, or Inf;
 * 3. all values are non-negative;
 * 4. at least min_integer_fraction of the values are close to integers;
 * 5. every sample has a library size > 0.
 */
// [[Rcpp::export]]
bool is_counts_matrix_cpp(SEXP initialized_matrix, bool verbose = false,
                          double integer_tol = 1e-8,
                          double min_integer_fraction = 0.95) {
  if (!std::isfinite(integer_tol) || integer_tol < 0.0) {
    stop("integer_tol must be a finite non-negative number.");
  }

  if (!std::isfinite(min_integer_fraction) || min_integer_fraction < 0.0 ||
      min_integer_fraction > 1.0) {
    stop("min_integer_fraction must be between 0 and 1.");
  }

  /*
   * Parse the beachmat external pointer.
   */
  Rtatami::BoundNumericPointer parsed(initialized_matrix);
  auto matrix = parsed->ptr;

  const std::size_t nr = static_cast<std::size_t>(matrix->nrow());

  const std::size_t nc = static_cast<std::size_t>(matrix->ncol());

  const std::size_t total_values = nr * nc;

  if (total_values == 0) {
    if (verbose) {
      cli_emit("info", "Result: FALSE; matrix is empty");
    }

    return false;
  }

  /*
   * Library size of each sample.
   */
  std::vector<double> library_sizes(nc, 0.0);

  std::size_t nonnegative_count = 0;
  std::size_t integer_like_count = 0;
  std::size_t zero_count = 0;

  /*
   * Read column by column.
   *
   * For sparse matrices, dense_column() returns unstored elements as 0;
   * for a DelayedMatrix it requests the data of the corresponding column.
   */
  auto accessor = matrix->dense_column();
  std::vector<double> buffer(nr);

  for (std::size_t col = 0; col < nc; ++col) {
    auto values = accessor->fetch(col, buffer.data());

    double sample_sum = 0.0;

    for (std::size_t row = 0; row < nr; ++row) {
      const double value = values[row];

      /*
       * std::isfinite() excludes all of:
       * - NA_real_
       * - NaN
       * - +Inf
       * - -Inf
       */
      if (!std::isfinite(value)) {
        if (verbose) {
          cli_emit("info", "Result: FALSE; matrix contains NA, NaN, or Inf");
        }

        return false;
      }

      if (value >= 0.0) {
        ++nonnegative_count;
      }

      if (value == 0.0) {
        ++zero_count;
      }

      /*
       * Use fabs(value - round(value)) to test whether the value is close to
       * an integer.
       */
      const double rounded = std::round(value);

      if (std::fabs(value - rounded) <= integer_tol) {
        ++integer_like_count;
      }

      /*
       * Keep the semantics of the original implementation:
       * the library size is the sum of all values in each column.
       */
      sample_sum += value;
    }

    library_sizes[col] = sample_sum;
  }

  /*
   * Check non-negativity.
   */
  const double nonnegative_fraction = static_cast<double>(nonnegative_count) /
                                      static_cast<double>(total_values);

  if (nonnegative_fraction < 1.0) {
    if (verbose) {
      cli_emit("info", "Result: FALSE; matrix contains negative values");
    }

    return false;
  }

  /*
   * Check the fraction of integer-like values.
   */
  const double integer_fraction = static_cast<double>(integer_like_count) /
                                  static_cast<double>(total_values);

  if (integer_fraction < min_integer_fraction) {
    if (verbose) {
      cli_emit("info", "Result: FALSE; insufficient integer-like values "
                       "(integer fraction " +
                           std::to_string(integer_fraction) + " < threshold " +
                           std::to_string(min_integer_fraction) + ")");
    }

    return false;
  }

  /*
   * Check the library size of each sample.
   *
   * Since Inf/NaN have already been excluded above, this mainly checks for:
   * - total <= 0;
   * - non-finite results, which should not occur in theory.
   */
  for (std::size_t col = 0; col < nc; ++col) {
    const double library_size = library_sizes[col];

    if (!std::isfinite(library_size) || library_size <= 0.0) {
      if (verbose) {
        cli_emit("info", "Result: FALSE; at least one sample has a "
                         "library size <= 0");
      }

      return false;
    }
  }

  /*
   * Success.
   */
  if (verbose) {
    const double zero_fraction =
        static_cast<double>(zero_count) / static_cast<double>(total_values);

    cli_emit("success",
             "Result: TRUE; matrix is likely a counts matrix "
             "(non-negative fraction " +
                 std::to_string(nonnegative_fraction) + ", integer fraction " +
                 std::to_string(integer_fraction) + ", zero fraction " +
                 std::to_string(zero_fraction) + ")");
  }

  return true;
}