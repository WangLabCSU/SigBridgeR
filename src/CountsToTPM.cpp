#include "Rtatami.h"
#include <Rcpp.h>

#include <cfloat>
#include <cmath>
#include <cstddef>

// inv_length: 每行基因的 1 / gene_length，由 R 端用 match() 对齐后传入。
// [[Rcpp::export]]
Rcpp::NumericMatrix counts_to_tpm_cpp(SEXP initialized_counts,
                                      const Rcpp::NumericVector &inv_length) {
  Rtatami::BoundNumericPointer parsed(initialized_counts);
  const auto &matrix = parsed->ptr;

  const int nrow = matrix->nrow();
  const int ncol = matrix->ncol();

  if (nrow <= 0 || ncol <= 0) {
    Rcpp::stop("counts must have positive nrow and ncol.");
  }
  if (inv_length.size() != static_cast<R_xlen_t>(nrow)) {
    Rcpp::stop("length(inv_length) must equal nrow(counts).");
  }

  // 所有元素都会被写入，无需清零
  Rcpp::NumericMatrix output = Rcpp::no_init(nrow, ncol);
  double *const out = REAL(output);
  const double *const inv = REAL(inv_length);
  const std::size_t n = static_cast<std::size_t>(nrow);

  auto accessor = matrix->dense_column();

  for (int col = 0; col < ncol; ++col) {
    if ((col & 63) == 0) {
      Rcpp::checkUserInterrupt();
    }

    double *const dst = out + static_cast<std::size_t>(col) * n;

    // 让 fetch 直接以输出列作缓冲区；若矩阵返回内部指针，则 src != dst
    const double *const src = accessor->fetch(col, dst);

    // 单次遍历：校验 + RPK + 分母（无分支，便于编译器优化）
    double denom = 0.0;
    bool bad = false;
    for (std::size_t i = 0; i < n; ++i) {
      const double v = src[i];
      bad |= !(v >= 0.0 && v <= DBL_MAX); // NaN / Inf / 负数都会失败
      const double rpk = v * inv[i];
      dst[i] = rpk;
      denom += rpk;
    }

    if (bad) {
      Rcpp::stop("counts must be finite and non-negative (column %d).",
                 col + 1);
    }
    if (!std::isfinite(denom)) {
      Rcpp::stop("invalid TPM denominator (column %d).", col + 1);
    }

    // denom == 0 说明整列为 0，dst 已全为 0，无需处理
    if (denom > 0.0) {
      const double factor = 1e6 / denom;
      for (std::size_t i = 0; i < n; ++i) {
        dst[i] *= factor;
      }
    }
  }

  return output;
}