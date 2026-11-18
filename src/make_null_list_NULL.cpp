#include <Rcpp.h>
using namespace Rcpp;

static SEXP make_null_list_NULL_rec(SEXP x) {
  // not a list, return as-is
  if (TYPEOF(x) != VECSXP) {
    return x;
  }
  // empty list() becomes NULL
  if (XLENGTH(x) == 0) {
    return R_NilValue;
  }
  R_xlen_t n = XLENGTH(x);
  List out(n);
  // preserve names
  SEXP names = Rf_getAttrib(x, R_NamesSymbol);
  if (names != R_NilValue) {
    out.attr("names") = names;
  }
  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = make_null_list_NULL_rec(VECTOR_ELT(x, i));
  }
  return out;
}

//' Replace Empty Lists with NULL Recursively
//'
//' @param x Any R object; lists are processed recursively.
//'
//' @return `x` with every empty list replaced by `NULL`; non-list objects are
//'   returned unchanged and names of non-empty lists are kept.
//'
//' @keywords internal
// [[Rcpp::export]]
SEXP make_null_list_NULL(SEXP x) { return make_null_list_NULL_rec(x); }
