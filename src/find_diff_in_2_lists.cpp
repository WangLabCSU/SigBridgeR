// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

using namespace Rcpp;

static const int IDENTICAL_FLAGS = 47;

static inline bool identical_cpp(SEXP x, SEXP y)
{
    return R_compute_identical(x, y, IDENTICAL_FLAGS);
}

static inline bool is_list_cpp(SEXP x)
{
    return TYPEOF(x) == VECSXP;
}

static inline SEXP no_diff()
{
    return Rf_ScalarLogical(1);
}

static inline bool is_no_diff(SEXP x)
{
    return TYPEOF(x) == LGLSXP &&
           XLENGTH(x) == 1 &&
           LOGICAL(x)[0] == 1;
}

static inline SEXP diff_pair(SEXP x, SEXP y, SEXP name)
{
    return List::create(
               _["name"] = name,
               _["x"] = x,
               _["y"] = y
           );
}

static inline bool has_names(SEXP x)
{
    SEXP nms = Rf_getAttrib(x, R_NamesSymbol);
    return nms != R_NilValue &&
           TYPEOF(nms) == STRSXP &&
           XLENGTH(nms) == XLENGTH(x);
}

static inline std::string get_name_str(SEXP nms, R_xlen_t i)
{
    if (nms == R_NilValue) return "";
    SEXP s = STRING_ELT(nms, i);
    if (s == NA_STRING) return "";
    return std::string(CHAR(s));
}

static inline SEXP get_name_sexp(SEXP nms, R_xlen_t i)
{
    if (nms == R_NilValue) return R_NilValue;
    SEXP s = STRING_ELT(nms, i);
    if (s == NA_STRING || std::strlen(CHAR(s)) == 0)
    {
        return R_NilValue;
    }
    return Rf_mkString(CHAR(s));
}

static SEXP find_diff_rec(SEXP x, SEXP y, SEXP current_name)
{
    if (x == y)
    {
        return no_diff();
    }
    if (is_list_cpp(x) && is_list_cpp(y))
    {
        R_xlen_t nx = XLENGTH(x);
        R_xlen_t ny = XLENGTH(y);
        bool use_names = has_names(x) && has_names(y);
        if (use_names)
        {
            SEXP xnms = Rf_getAttrib(x, R_NamesSymbol);
            SEXP ynms = Rf_getAttrib(y, R_NamesSymbol);
            std::unordered_map<std::string, R_xlen_t> ypos;
            for (R_xlen_t j = 0; j < ny; ++j)
            {
                std::string nm = get_name_str(ynms, j);
                if (!nm.empty() && ypos.find(nm) == ypos.end())
                {
                    ypos[nm] = j;
                }
            }
            std::unordered_set<std::string> seen_x;
            for (R_xlen_t i = 0; i < nx; ++i)
            {
                std::string nm = get_name_str(xnms, i);
                SEXP nm_sexp = get_name_sexp(xnms, i);
                if (nm.empty())
                {
                    if (i >= ny)
                    {
                        return diff_pair(VECTOR_ELT(x, i), R_NilValue, current_name);
                    }
                    SEXP ans = find_diff_rec(
                                   VECTOR_ELT(x, i),
                                   VECTOR_ELT(y, i),
                                   current_name
                               );
                    if (!is_no_diff(ans)) return ans;
                    continue;
                }
                seen_x.insert(nm);
                auto it = ypos.find(nm);
                if (it == ypos.end())
                {
                    return diff_pair(VECTOR_ELT(x, i), R_NilValue, nm_sexp);
                }
                SEXP ans = find_diff_rec(
                               VECTOR_ELT(x, i),
                               VECTOR_ELT(y, it->second),
                               nm_sexp
                           );
                if (!is_no_diff(ans)) return ans;
            }
            for (R_xlen_t j = 0; j < ny; ++j)
            {
                std::string nm = get_name_str(ynms, j);
                SEXP nm_sexp = get_name_sexp(ynms, j);
                if (!nm.empty() && seen_x.find(nm) == seen_x.end())
                {
                    return diff_pair(R_NilValue, VECTOR_ELT(y, j), nm_sexp);
                }
                if (nm.empty() && j >= nx)
                {
                    return diff_pair(R_NilValue, VECTOR_ELT(y, j), current_name);
                }
            }
            return no_diff();
        }
        R_xlen_t n = std::min(nx, ny);
        for (R_xlen_t i = 0; i < n; ++i)
        {
            SEXP ans = find_diff_rec(
                           VECTOR_ELT(x, i),
                           VECTOR_ELT(y, i),
                           current_name
                       );
            if (!is_no_diff(ans)) return ans;
        }
        if (nx > ny)
        {
            return diff_pair(VECTOR_ELT(x, n), R_NilValue, current_name);
        }
        if (ny > nx)
        {
            return diff_pair(R_NilValue, VECTOR_ELT(y, n), current_name);
        }
        return no_diff();
    }
    if (identical_cpp(x, y))
    {
        return no_diff();
    }
    return diff_pair(x, y, current_name);
}

// [[Rcpp::export]]
SEXP find_diff_in_2_lists(SEXP x, SEXP y)
{
    if (!is_list_cpp(x) || !is_list_cpp(y))
    {
        stop("x and y must both be lists.");
    }
    return find_diff_rec(x, y, R_NilValue);
}