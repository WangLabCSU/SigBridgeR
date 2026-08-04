#include <Rcpp.h>
using namespace Rcpp;

static SEXP get_element(SEXP x, int i)
{
    int idx = i - 1;
    if (TYPEOF(x) == VECSXP)
    {
        return VECTOR_ELT(x, idx);
    }
    switch (TYPEOF(x))
    {
        case INTSXP:
            return Rf_ScalarInteger(INTEGER(x)[idx]);
        case REALSXP:
            return Rf_ScalarReal(REAL(x)[idx]);
        case LGLSXP:
            return Rf_ScalarLogical(LOGICAL(x)[idx]);
        case STRSXP:
            return Rf_ScalarString(STRING_ELT(x, idx));
        default:
            stop("Unsupported rhs type");
    }
    return R_NilValue;
}

static SEXP get_path(SEXP rhs, IntegerVector path)
{
    SEXP cur = rhs;
    for (int i = 0; i < path.size(); ++i)
    {
        cur = get_element(cur, path[i]);
    }
    return cur;
}

// [[Rcpp::export]]
void fast_assign_plan(List plan, SEXP rhs, Environment env)
{
    for (int i = 0; i < plan.size(); ++i)
    {
        List item = plan[i];
        std::string type = as<std::string>(item["type"]);
        if (type == "ignore") continue;
        std::string name = as<std::string>(item["name"]);
        IntegerVector path = item["path"];
        SEXP value = get_path(rhs, path);
        env.assign(name, value);
    }
}