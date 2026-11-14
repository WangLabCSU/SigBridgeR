// src/resolve_param_case.cpp
#include "SigBridgeR.h"
#include <Rcpp.h>

#include <cctype>
#include <string>
#include <vector>

namespace {

inline std::string ascii_lower(const char *s) {
  std::string out(s);
  for (char &c : out) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }
  return out;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List resolve_param_case(const Rcpp::List &dots,
                              const Rcpp::CharacterVector &canonical,
                              Rcpp::Environment env, bool lower_inner = false) {
  SEXP nms = Rf_getAttrib(dots, R_NamesSymbol);
  if (Rf_isNull(nms)) {
    return dots; // 快速路径：没有任何具名参数
  }

  const R_xlen_t n = dots.size();
  const R_xlen_t m = canonical.size();

  std::vector<std::string> canon_lower;
  canon_lower.reserve(m);
  for (R_xlen_t j = 0; j < m; ++j) {
    canon_lower.push_back(ascii_lower(CHAR(STRING_ELT(canonical, j))));
  }

  std::vector<char> is_hit(n, 0);
  R_xlen_t n_hit = 0;

  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP nm = STRING_ELT(nms, i);
    if (nm == NA_STRING || CHAR(nm)[0] == '\0') {
      continue;
    }

    const std::string low = ascii_lower(CHAR(nm));
    R_xlen_t k = -1;
    for (R_xlen_t j = 0; j < m; ++j) {
      if (canon_lower[j] == low) {
        k = j;
        break;
      }
    }
    if (k < 0) {
      continue;
    }

    SEXP val_sexp = VECTOR_ELT(dots, i);
    if (TYPEOF(val_sexp) != VECSXP) {
      Rcpp::stop("参数 `%s` 必须是 list，当前为 %s", CHAR(nm),
                 Rf_type2char(TYPEOF(val_sexp)));
    }
    Rcpp::List val(val_sexp);

    if (lower_inner) {
      SEXP vn = Rf_getAttrib(val_sexp, R_NamesSymbol);
      if (!Rf_isNull(vn)) {
        // 只复制外层向量，元素本身共享，不做深拷贝
        val = Rcpp::List(Rf_shallow_duplicate(val_sexp));
        const R_xlen_t len = Rf_xlength(vn);
        Rcpp::CharacterVector new_names(len);
        for (R_xlen_t t = 0; t < len; ++t) {
          SEXP s = STRING_ELT(vn, t);
          new_names[t] = (s == NA_STRING)
                             ? NA_STRING
                             : Rf_mkChar(ascii_lower(CHAR(s)).c_str());
        }
        val.attr("names") = new_names;
      }
    }

    const char *target = CHAR(STRING_ELT(canonical, k));

    // 直接按引用取当前值（会强制 promise），不经过 R 层 get()
    SEXP cur_sexp = env.get(target);
    Rcpp::List cur;
    if (!Rf_isNull(cur_sexp)) {
      if (TYPEOF(cur_sexp) != VECSXP) {
        Rcpp::stop("参数 `%s` 必须是 list", target);
      }
      cur = Rcpp::List(cur_sexp);
    }

    // 以规范名参数为准：等价于 modifyList(val, cur)
    env.assign(target, sigbridger::modify_list(val, cur));

    is_hit[i] = 1;
    ++n_hit;
  }

  if (n_hit == 0) {
    return dots;
  }

  Rcpp::List out(n - n_hit);
  Rcpp::CharacterVector out_names(n - n_hit);
  R_xlen_t w = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    if (!is_hit[i]) {
      out[w] = VECTOR_ELT(dots, i);
      out_names[w] = STRING_ELT(nms, i);
      ++w;
    }
  }
  out.attr("names") = out_names;
  return out;
}