#include "Rtatami.h"
#include <Rcpp.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

using namespace Rcpp;

namespace {
constexpr int kOk = 0;
constexpr int kNaNotAllowed = 1; // 遇到 NA 且 na_as_zero = FALSE
constexpr int kNegative = 2;
constexpr int kInfinite = 3;
constexpr int kTotalOverflow = 4; // 列总和溢出为 Inf
constexpr int kTotalTooSmall = 5; // 1e6 / total 溢出
} // namespace

// 返回 list(tpm, status, row, col, na_count, zero_cols)
// status != 0 时 tpm 为 NULL，row/col 为首个问题位置（1-based，col 可能无 row）
// [[Rcpp::export]]
List fpkm_to_tpm_cpp(SEXP initialized_fpkm, bool na_as_zero = true) {
  Rtatami::BoundNumericPointer parsed(initialized_fpkm);
  const auto &matrix = parsed->ptr;
  const int nr = matrix->nrow();
  const int nc = matrix->ncol();

  NumericMatrix output(no_init(nr, nc));
  double *const out_base = output.begin();

  int status = kOk;
  int bad_row = 0;
  int bad_col = 0;
  double na_count = 0.0;      // 用 double 计数，避免超过 2^31 个 NA 时溢出
  std::vector<int> zero_cols; // 总和为 0 的列（1-based）

  if (nr > 0 && nc > 0) {
    constexpr double kMax = std::numeric_limits<double>::max();
    auto ext = matrix->dense_column();
    std::vector<double> buffer(static_cast<std::size_t>(nr));

    // 单趟：读一次 -> 校验 -> 写入输出列 -> 原地缩放（该列仍在缓存中）
    for (int c = 0; c < nc; ++c) {
      if ((c & 0xFF) == 0) {
        checkUserInterrupt();
      }

      const double *values = ext->fetch(c, buffer.data());
      double *const out = out_base + static_cast<R_xlen_t>(c) * nr;
      double total = 0.0;

      for (int r = 0; r < nr; ++r) {
        const double v = values[r];
        if (!(v >= 0.0 && v <= kMax)) { // NaN / 负值 / ±Inf，一次比较搞定热路径
          if (std::isnan(v)) {
            if (!na_as_zero) {
              status = kNaNotAllowed;
              bad_row = r + 1;
              bad_col = c + 1;
              break;
            }
            ++na_count;
            out[r] = 0.0;
            continue;
          }
          status = std::isinf(v) ? kInfinite : kNegative;
          bad_row = r + 1;
          bad_col = c + 1;
          break;
        }
        out[r] = v;
        total += v;
      }
      if (status != kOk) {
        break;
      }

      if (!std::isfinite(total)) {
        status = kTotalOverflow;
        bad_col = c + 1;
        break;
      }
      if (total > 0.0) {
        const double scale = 1e6 / total;
        if (!std::isfinite(scale)) {
          status = kTotalTooSmall;
          bad_col = c + 1;
          break;
        }
        for (int r = 0; r < nr; ++r) {
          out[r] *= scale;
        }
      } else {
        // 全部非负且总和为 0：整列已经是 0（NA 也已写成 0）
        zero_cols.push_back(c + 1);
      }
    }
  }

  return List::create(
      Named("tpm") = (status == kOk ? static_cast<SEXP>(output) : R_NilValue),
      Named("status") = status, Named("row") = bad_row, Named("col") = bad_col,
      Named("na_count") = na_count,
      Named("zero_cols") = IntegerVector(zero_cols.begin(), zero_cols.end()));
}