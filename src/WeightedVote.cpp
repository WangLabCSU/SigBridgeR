#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

namespace {

struct TieRandom {
  static inline bool positive_wins() { return R::unif_rand() < 0.5; }
};
struct TieFirst {
  static inline bool positive_wins() { return true; }
};
struct TieLast {
  static inline bool positive_wins() { return false; }
};

constexpr double kRelTol = 1e-12;

template <typename TiePolicy>
CharacterVector weighted_vote_impl(const CharacterMatrix &vote_data,
                                   const NumericVector &weights) {
  const R_xlen_t nr = vote_data.nrow();
  const R_xlen_t nc = vote_data.ncol();
  if (weights.size() != nc) {
    stop("Length of `weights` must equal ncol(vote_data).");
  }

  // 1) 一次性校验并压缩有效列（NA/NaN/0 不贡献）
  std::vector<R_xlen_t> cols;
  std::vector<double> w;
  cols.reserve(nc);
  w.reserve(nc);
  for (R_xlen_t j = 0; j < nc; ++j) {
    const double wj = weights[j];
    if (ISNAN(wj) || wj == 0.0)
      continue;
    if (!R_FINITE(wj) || wj < 0.0) {
      stop("`weights` must be finite and non-negative (NA/NaN allowed).");
    }
    cols.push_back(j);
    w.push_back(wj);
  }

  // 2) 预先创建并保护结果用的 CHARSXP，后续用指针比较
  Rcpp::Shield<SEXP> pos_sx(Rf_mkChar("Positive"));
  Rcpp::Shield<SEXP> oth_sx(Rf_mkChar("Other"));
  SEXP kPos = pos_sx;
  SEXP kOth = oth_sx;

  // 3) 列优先累加，顺序访问内存
  SEXP vs = vote_data;
  const SEXP *v = STRING_PTR_RO(vs);
  std::vector<double> sp(nr, 0.0), so(nr, 0.0);
  for (std::size_t k = 0; k < cols.size(); ++k) {
    const SEXP *col = v + static_cast<R_xlen_t>(cols[k]) * nr;
    const double wk = w[k];
    for (R_xlen_t i = 0; i < nr; ++i) {
      SEXP s = col[i];
      if (s == NA_STRING)
        continue;
      if (s == kPos)
        sp[i] += wk;
      else
        so[i] += wk;
    }
  }

  // 4) 按行序决定胜负，保证随机平局的 RNG 消耗顺序与逐行版本一致
  CharacterVector out(nr);
  for (R_xlen_t i = 0; i < nr; ++i) {
    const double p = sp[i], o = so[i];
    bool pos_wins;
    if (std::fabs(p - o) <= kRelTol * std::max(p, o)) {
      pos_wins = TiePolicy::positive_wins();
    } else {
      pos_wins = p > o;
    }
    SET_STRING_ELT(out, i, pos_wins ? kPos : kOth);
  }
  return out;
}

} // namespace

// [[Rcpp::export(rng = false)]]
CharacterVector weighted_vote_cpp(CharacterMatrix vote_data,
                                  NumericVector weights, int ties_method) {
  switch (ties_method) {
  case 0: {
    Rcpp::RNGScope scope; // 仅随机平局需要
    return weighted_vote_impl<TieRandom>(vote_data, weights);
  }
  case 1:
    return weighted_vote_impl<TieFirst>(vote_data, weights);
  case 2:
    return weighted_vote_impl<TieLast>(vote_data, weights);
  default:
    stop("Invalid `ties_method`.");
  }
}