#include "Rtatami.h"
#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace Rcpp;

enum class AggMethod { Max, Sum, Mean, Median, First };

static AggMethod parse_method(const std::string &m) {
  if (m == "max")
    return AggMethod::Max;
  if (m == "sum")
    return AggMethod::Sum;
  if (m == "mean")
    return AggMethod::Mean;
  if (m == "median")
    return AggMethod::Median;
  if (m == "first")
    return AggMethod::First;
  stop("Unsupported method: " + m);
  return AggMethod::Max;
}

// CSR 布局：组 g 的成员为 member[offset[g] .. offset[g+1])，组内保持原始升序
struct Groups {
  std::vector<int> offset;
  std::vector<int> member;
};

static Groups build_groups(const IntegerVector &gid, int n, int ng) {
  if (gid.size() != n) {
    stop("Length of group ids must match the dimension being aggregated.");
  }
  Groups g;
  g.offset.assign(static_cast<size_t>(ng) + 1, 0);

  for (int i = 0; i < n; ++i) {
    const int k = gid[i];
    if (k == NA_INTEGER || k < 1 || k > ng) {
      stop("Invalid group id.");
    }
    ++g.offset[k]; // 先把 "组 k-1 的个数" 记在 k 位置
  }
  for (int k = 0; k < ng; ++k) {
    g.offset[k + 1] += g.offset[k];
  }

  std::vector<int> cursor(g.offset.begin(), g.offset.end() - 1);
  g.member.resize(n);
  for (int i = 0; i < n; ++i) {
    g.member[cursor[gid[i] - 1]++] = i;
  }
  return g;
}

static inline double median_inplace(std::vector<double> &v) {
  if (v.empty())
    return NA_REAL;
  const size_t n = v.size(), mid = n / 2;
  std::nth_element(v.begin(), v.begin() + mid, v.end());
  const double hi = v[mid];
  if (n & 1)
    return hi;
  const double lo = *std::max_element(v.begin(), v.begin() + mid);
  return (lo + hi) / 2.0;
}

// 对 v[idx[0..k)] 聚合（rows 版使用，tmp 供 median 复用）
static inline double aggregate_gather(const double *v, const int *idx, int k,
                                      AggMethod m, std::vector<double> &tmp) {
  switch (m) {
  case AggMethod::First:
    return v[idx[0]];
  case AggMethod::Sum: {
    double s = 0.0;
    for (int t = 0; t < k; ++t) {
      const double x = v[idx[t]];
      if (!std::isnan(x))
        s += x;
    }
    return s;
  }
  case AggMethod::Mean: {
    double s = 0.0;
    int cnt = 0;
    for (int t = 0; t < k; ++t) {
      const double x = v[idx[t]];
      if (!std::isnan(x)) {
        s += x;
        ++cnt;
      }
    }
    return cnt == 0 ? R_NaN : s / cnt;
  }
  case AggMethod::Max: {
    double best = NA_REAL;
    for (int t = 0; t < k; ++t) {
      const double x = v[idx[t]];
      if (!std::isnan(x) && (std::isnan(best) || x > best))
        best = x;
    }
    return best;
  }
  case AggMethod::Median: {
    tmp.clear();
    for (int t = 0; t < k; ++t) {
      const double x = v[idx[t]];
      if (!std::isnan(x))
        tmp.push_back(x);
    }
    return median_inplace(tmp);
  }
  }
  return NA_REAL;
}

// [[Rcpp::export]]
NumericMatrix aggregate_dup_cols_cpp(SEXP initialized_matrix,
                                     IntegerVector group_id, int n_groups,
                                     std::string method) {
  Rtatami::BoundNumericPointer parsed(initialized_matrix);
  const auto &mat = parsed->ptr;
  const int nr = mat->nrow();
  const int nc = mat->ncol();

  const AggMethod agg = parse_method(method);
  const Groups grp = build_groups(group_id, nc, n_groups);

  NumericMatrix output(Rcpp::no_init(nr, n_groups));
  double *out = REAL(output);

  auto ext = mat->dense_column();
  std::vector<double> buf(nr);
  std::vector<int> cnt;           // mean 用
  std::vector<double> block, tmp; // median 用

  for (int g = 0; g < n_groups; ++g) {
    if ((g & 0xFF) == 0)
      Rcpp::checkUserInterrupt();

    double *dst = out + static_cast<size_t>(g) * nr;
    const int *mb = grp.member.data() + grp.offset[g];
    const int k = grp.offset[g + 1] - grp.offset[g];

    // 单例 / first：原样拷贝，行为与 rows 版一致
    if (k == 1 || agg == AggMethod::First) {
      const double *p = ext->fetch(mb[0], buf.data());
      std::copy(p, p + nr, dst);
      continue;
    }

    switch (agg) {
    case AggMethod::Sum: {
      std::fill(dst, dst + nr, 0.0);
      for (int t = 0; t < k; ++t) {
        const double *p = ext->fetch(mb[t], buf.data());
        for (int i = 0; i < nr; ++i)
          if (!std::isnan(p[i]))
            dst[i] += p[i];
      }
      break;
    }
    case AggMethod::Mean: {
      std::fill(dst, dst + nr, 0.0);
      cnt.assign(nr, 0);
      for (int t = 0; t < k; ++t) {
        const double *p = ext->fetch(mb[t], buf.data());
        for (int i = 0; i < nr; ++i)
          if (!std::isnan(p[i])) {
            dst[i] += p[i];
            ++cnt[i];
          }
      }
      for (int i = 0; i < nr; ++i)
        dst[i] = cnt[i] ? dst[i] / cnt[i] : R_NaN;
      break;
    }
    case AggMethod::Max: {
      std::fill(dst, dst + nr, NA_REAL);
      for (int t = 0; t < k; ++t) {
        const double *p = ext->fetch(mb[t], buf.data());
        for (int i = 0; i < nr; ++i)
          if (!std::isnan(p[i]) && (std::isnan(dst[i]) || p[i] > dst[i]))
            dst[i] = p[i];
      }
      break;
    }
    case AggMethod::Median: {
      block.resize(static_cast<size_t>(k) * nr);
      for (int t = 0; t < k; ++t) {
        const double *p = ext->fetch(mb[t], buf.data());
        std::copy(p, p + nr, block.begin() + static_cast<size_t>(t) * nr);
      }
      for (int i = 0; i < nr; ++i) {
        tmp.clear();
        for (int t = 0; t < k; ++t) {
          const double x = block[static_cast<size_t>(t) * nr + i];
          if (!std::isnan(x))
            tmp.push_back(x);
        }
        dst[i] = median_inplace(tmp);
      }
      break;
    }
    case AggMethod::First:
      break; // 已在上面处理
    }
  }
  return output;
}

// [[Rcpp::export]]
NumericMatrix aggregate_dup_rows_cpp(SEXP initialized_matrix,
                                     IntegerVector group_id, int n_groups,
                                     std::string method) {
  Rtatami::BoundNumericPointer parsed(initialized_matrix);
  const auto &mat = parsed->ptr;
  const int nr = mat->nrow();
  const int nc = mat->ncol();

  const AggMethod agg = parse_method(method);
  const Groups grp = build_groups(group_id, nr, n_groups);

  NumericMatrix output(Rcpp::no_init(n_groups, nc));
  double *out = REAL(output);

  auto ext = mat->dense_column();
  std::vector<double> buf(nr), tmp;

  for (int j = 0; j < nc; ++j) {
    if ((j & 0xFF) == 0)
      Rcpp::checkUserInterrupt();

    const double *v = ext->fetch(j, buf.data());
    double *dst = out + static_cast<size_t>(j) * n_groups;

    for (int g = 0; g < n_groups; ++g) {
      const int *mb = grp.member.data() + grp.offset[g];
      const int k = grp.offset[g + 1] - grp.offset[g];
      dst[g] = (k == 1) ? v[mb[0]] : aggregate_gather(v, mb, k, agg, tmp);
    }
  }
  return output;
}