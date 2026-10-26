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

  Rtatami::BoundNumericPointer parsed(initialized_matrix);
  const auto &matrix = parsed->ptr;
  const auto nr = matrix->nrow();
  const auto nc = matrix->ncol();

  auto fail = [&](const char *msg) {
    if (verbose)
      cli_emit("info", std::string(msg));
    return false;
  };

  if (nr == 0 || nc == 0)
    return fail("Result: FALSE; matrix is empty");

  std::size_t nonzero = 0;      // 非零元素数
  std::size_t integer_like = 0; // 非零元素中接近整数的个数

  // 返回 0 = 正常；1 = 非有限；2 = 负值。零值最先短路。
  auto classify = [&](double v) -> int {
    if (v == 0.0)
      return 0;
    if (!std::isfinite(v))
      return 1;
    if (v < 0.0)
      return 2;
    ++nonzero;
    if (std::fabs(v - std::nearbyint(v)) <= integer_tol)
      ++integer_like;
    return 0;
  };

  auto reject = [&](int code) {
    return code == 1 ? fail("Result: FALSE; matrix contains NA, NaN, or Inf")
                     : fail("Result: FALSE; matrix contains negative values");
  };

  // 值非负且有限时，library size > 0 等价于该列至少有一个正值
  std::size_t empty_cols = 0;

  if (matrix->is_sparse()) {
    tatami::Options opt;
    opt.sparse_extract_index = false; // 只需要值
    auto ext = matrix->sparse_column(opt);
    std::vector<double> vbuf(static_cast<std::size_t>(nr));

    for (int c = 0; c < nc; ++c) {
      if ((c & 0xFF) == 0)
        Rcpp::checkUserInterrupt();
      auto range = ext->fetch(c, vbuf.data(), nullptr);
      std::size_t before = nonzero;
      for (int k = 0; k < range.number; ++k) {
        if (int code = classify(range.value[k]))
          return reject(code);
      }
      if (nonzero == before)
        ++empty_cols;
    }
  } else {
    auto ext = matrix->dense_column();
    std::vector<double> buf(static_cast<std::size_t>(nr));

    for (int c = 0; c < nc; ++c) {
      if ((c & 0xFF) == 0)
        Rcpp::checkUserInterrupt();
      const double *values = ext->fetch(c, buf.data());
      std::size_t before = nonzero;
      for (int r = 0; r < nr; ++r) {
        if (int code = classify(values[r]))
          return reject(code);
      }
      if (nonzero == before)
        ++empty_cols;
    }
  }

  if (empty_cols > 0) {
    return fail("Result: FALSE; at least one sample has a library size <= 0");
  }

  // nonzero > 0 在这里必然成立（每列至少有一个正值）
  const double integer_fraction =
      static_cast<double>(integer_like) / static_cast<double>(nonzero);
  if (integer_fraction < min_integer_fraction) {
    if (verbose) {
      cli_emit("info", "Result: FALSE; insufficient integer-like values among "
                       "non-zero entries (integer fraction " +
                           std::to_string(integer_fraction) + " < threshold " +
                           std::to_string(min_integer_fraction) + ")");
    }
    return false;
  }

  if (verbose) {
    const double density = static_cast<double>(nonzero) /
                           (static_cast<double>(nr) * static_cast<double>(nc));
    cli_emit("success", "Result: TRUE; matrix is likely a counts matrix "
                        "(integer fraction among non-zero " +
                            std::to_string(integer_fraction) + ", density " +
                            std::to_string(density) + ")");
  }
  return true;
}