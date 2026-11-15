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



//' Resolve `...` Names Against Canonical Parameter Names
//'
//' Matches the names of a captured `...` list against a set of canonical
//' parameter names, ignoring case, and writes the matched values into the
//' calling frame. Matched elements are dropped from the returned list, so the
//' caller can forward what is left as genuine `...`.
//'
//' @details
//' This is a helper for functions that accept a fixed set of named
//' sub-argument lists (for example `bench_config`) but want to tolerate
//' `BenchRounds` as well as `bench_rounds`.
//'
//' * Matching is case-insensitive on the **ASCII** range only, and is
//'   otherwise exact: `foo_bar` does not match `foo-bar`. If `canonical`
//'   contains several names that differ only in case, the first one wins.
//' * Each matched value must be a plain list, otherwise an error is raised
//'   naming the offending argument.
//' * The current value of the canonical parameter is read from `env` (which
//'   forces its promise) and merged with the supplied value as
//'   `modify_list(supplied, current)`: fields present in both take the value
//'   from `env`, fields only in the supplied list are added, nested lists are
//'   merged recursively, and `NULL` values are kept. The result is then
//'   assigned back to the canonical name in `env`.
//' * When several supplied names match the same canonical parameter, they are
//'   processed in order, so the last one wins.
//' * Elements of `dots` that are unnamed, or whose name is `NA` or `""`, are
//'   never matched and are returned unchanged.
//' * If `dots` has no names at all, it is returned as is without any work.
//' * Unmatched elements and untouched values are shared with the input, not
//'   copied.
//'
//' @param dots A list, typically `list(...)` from the calling function. Its
//'   names are compared against `canonical`.
//' @param canonical A character vector of canonical parameter names, typically
//'   `rlang::fn_fmls_names()`. The exact spelling of a matched entry is used as
//'   the name of the variable assigned in `env`.
//' @param env An environment, typically `rlang::current_env()`, in which the
//'   canonical parameters live. Matched parameters are written here; they are
//'   not part of the returned list.
//' @param lower_inner A logical flag. When `TRUE`, the inner names of every
//'   matched list are lowercased as well, so `list(Alpha = list(MaxIter = 1))`
//'   becomes `list(Alpha = list(maxiter = 1))`.
//'
//' @return A named list holding the elements of `dots` that were not matched,
//'   in their original order. Its length may be zero.
//'
//' @examples
//' \dontrun{
//' f <- function(alpha = list(), ...) {
//'   dots <- resolve_param_case(
//'     dots = list(...),
//'     canonical = rlang::fn_fmls_names(),
//'     env = rlang::current_env(),
//'     lower_inner = TRUE
//'   )
//'   list(alpha = alpha, dots = dots)
//' }
//'
//' # `Alpha` matches `alpha`; inner names are lowercased, `extra` is left over
//' f(Alpha = list(MaxIter = 10), extra = 1) |> str()
//'
//' # non-list values are rejected
//' try(f(Alpha = 10))
//' }
//'
//' @seealso [modify_list()]
//' @keywords internal
//'
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