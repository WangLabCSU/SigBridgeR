// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

static SEXP make_null_list_NULL_rec(SEXP x)
{
    // not a list, return as-is
    if (TYPEOF(x) != VECSXP)
    {
        return x;
    }
    // empty list() becomes NULL
    if (XLENGTH(x) == 0)
    {
        return R_NilValue;
    }
    R_xlen_t n = XLENGTH(x);
    List out(n);
    // preserve names
    SEXP names = Rf_getAttrib(x, R_NamesSymbol);
    if (names != R_NilValue)
    {
        out.attr("names") = names;
    }
    for (R_xlen_t i = 0; i < n; ++i)
    {
        out[i] = make_null_list_NULL_rec(VECTOR_ELT(x, i));
    }
    return out;
}

// [[Rcpp::export]]
SEXP make_null_list_NULL(SEXP x)
{
    return make_null_list_NULL_rec(x);
}
