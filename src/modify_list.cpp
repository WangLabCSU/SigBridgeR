#include "SigBridgeR.h"

// [[Rcpp::plugins(cpp17)]]

using namespace Rcpp;

//' Faster modifyList()
//' 
//' @inherit utils::modifyList
//' @param keep_null If TRUE, NULL elements in val become NULL elements in x. Otherwise, the corresponding element, if present, is deleted from x.
//' @keywords internal
// [[Rcpp::export]]
List modify_list(List x, List val, bool keep_null = true) {
  return sigbridger::modify_list(x, val, keep_null);
}

//' Merge Matching Sub-lists of `val` Against `x`
//'
//' For every element of `val`, if `x` has an element with the same name and
//' both are lists, the two are merged recursively with `modify_list()`;
//' otherwise the element of `val` is used unchanged.
//'
//' Note that the result is built from `val` only: its length and names are
//' those of `val`, and elements of `x` whose names do not appear in `val`
//' are **not** carried over. This differs from [utils::modifyList()].
//'
//' @details
//' * Names are matched exactly (case-sensitive). If `x` has duplicated names,
//'   the first match is used.
//' * Elements of `val` that are unnamed, have an `NA` name or an empty name
//'   are passed through unchanged.
//' * Unmodified elements are shared with the inputs, not deep-copied.
//' * Only the `names` attribute of `val` is kept on the result.
//'
//' @param x A list. Provides the base values for matching sub-lists.
//' @param val A list. Provides the new values.
//'
//' @return A list with the same length and names as `val`.
//'
//' @examples
//' \dontrun{
//' x <- list(
//'   a = list(p = 1, q = 2),
//'   b = 1,
//'   c = list(z = 0)
//' )
//' val <- list(
//'   a = list(q = 20, r = 30), # both lists: merged recursively
//'   b = list(new = TRUE),     # `x$b` is not a list: replaced by `val$b`
//'   d = 4                     # not in `x`: kept as is
//' )
//'
//' res <- SigBridgeR:::modify_list2(x, val)
//' str(res)
//' # res$a: p = 1, q = 20, r = 30
//' # res$b: list(new = TRUE)
//' # res$d: 4
//' # `x$c` is not in `val`, so it is not in the result
//' names(res)
//'
//' # unnamed elements are passed through
//' SigBridgeR:::modify_list2(list(1), list(list(a = 1), 2))
//'
//' # non-list inputs are rejected
//' try(SigBridgeR:::modify_list2(1:3, list()))
//' }
//' @keywords internal
//'
// [[Rcpp::export]]
Rcpp::List modify_list2(SEXP x, SEXP val) {
  // 用 TYPEOF 而不是 Rcpp::List 形参，避免非 list 被 Rcpp 静默转换
  if (!sigbridger::is_plain_list(x)) {
    Rcpp::stop("`x` must be a list.");
  }
  if (!sigbridger::is_plain_list(val)) {
    Rcpp::stop("`val` must be a list.");
  }

  const R_xlen_t n = Rf_xlength(val);
  Rcpp::List out(n);

  SEXP val_nms = Rf_getAttrib(val, R_NamesSymbol);
  SEXP x_nms = Rf_getAttrib(x, R_NamesSymbol);

  // 任何一方没有名字或为空，就不可能发生匹配
  const bool can_match =
      n > 0 && !Rf_isNull(val_nms) && !Rf_isNull(x_nms) && Rf_xlength(x) > 0;

  // 惰性构建 x 的名字索引：只有真正遇到 list 类型的 val 元素时才建
  std::unordered_map<std::string, R_xlen_t> index;
  bool index_built = false;

  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP v = VECTOR_ELT(val, i);

    if (can_match && TYPEOF(v) == VECSXP) {
      SEXP key = STRING_ELT(val_nms, i);

      if (key != NA_STRING && CHAR(key)[0] != '\0') {
        if (!index_built) {
          const R_xlen_t xn = Rf_xlength(x_nms);
          index.reserve(static_cast<size_t>(xn));
          for (R_xlen_t j = 0; j < xn; ++j) {
            SEXP nm = STRING_ELT(x_nms, j);
            if (nm == NA_STRING || CHAR(nm)[0] == '\0') {
              continue;
            }
            // emplace 不会覆盖已有键，因此保留首次出现，与 x[[key]] 一致
            index.emplace(Rf_translateCharUTF8(nm), j);
          }
          index_built = true;
        }

        auto it = index.find(Rf_translateCharUTF8(key));
        if (it != index.end()) {
          SEXP cur = VECTOR_ELT(x, it->second);
          if (TYPEOF(cur) == VECSXP) {
            // Rcpp::List(SEXP) 只是句柄，不复制数据
            Rcpp::List merged =
                sigbridger::modify_list(Rcpp::List(cur), Rcpp::List(v));
            SET_VECTOR_ELT(out, i, merged);
            continue;
          }
        }
      }
    }

    SET_VECTOR_ELT(out, i, v); // 共享引用，不做深拷贝
  }

  if (!Rf_isNull(val_nms)) {
    Rf_setAttrib(out, R_NamesSymbol, val_nms);
  }

  return out;
}